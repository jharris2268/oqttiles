#include "visalingham_whyatt.hpp"
#include <cmath>
#include <algorithm>
#include <pybind11/pybind11.h>
#include <limits>
#include <pybind11/stl.h>


#include "geos_wrapper.hpp"

double triangle_area(
    const std::tuple<double,double,double>& p1,
    const std::tuple<double,double,double>& p2,
    const std::tuple<double,double,double>& p3) {
    
    
    return abs(
        std::get<0>(p1)*(std::get<1>(p2)-std::get<1>(p3)) +
        std::get<0>(p2)*(std::get<1>(p3)-std::get<1>(p1)) +
        std::get<0>(p3)*(std::get<1>(p1)-std::get<1>(p2))
    )/2.;
}

std::tuple<size_t,size_t,double,size_t> find_and_drop_min(std::vector<std::pair<size_t,double>>& real_areas, size_t curr_len) {
    if (curr_len <= 2) {
        return std::make_tuple(0, 0, std::numeric_limits<double>::infinity(), curr_len);
    }
    auto it = std::min_element(real_areas.begin()+1, real_areas.begin()+curr_len-1, [](const auto& l, const auto& r) { return l.second<r.second; });
    size_t idx = it - real_areas.begin();
    size_t curr=it->first;
    double min_area=it->second;
    
    if (idx < (curr_len-1)) {
        std::copy(real_areas.begin() + idx + 1, real_areas.begin()+curr_len, real_areas.begin()+idx);
    }
    
    return std::make_tuple(idx,curr,min_area,curr_len-1);
}
    
    


VWSimplify::VWSimplify(const std::vector<std::pair<double,double>>& coords_in) {
    
    if (coords_in.size()<2) {
        throw std::domain_error("??");
    }
    
    coords_.reserve(coords_in.size());
    ordered_weights_.reserve(coords_in.size());    
    
    std::vector<std::pair<size_t,double>> real_areas;
    
    double inf=std::numeric_limits<double>::infinity();
    
    for (size_t i=0; i < coords_in.size(); i++) {
        
        coords_.push_back(std::make_tuple(coords_in[i].first,coords_in[i].second,inf));
        ordered_weights_.push_back(i);
    }
    
    real_areas.push_back(std::make_pair(0,inf));
    for (size_t i=1; i < coords_.size()-1; i++) {
        double a = triangle_area(coords_[i-1],coords_[i],coords_[i+1]);
        std::get<2>(coords_[i]) = a;
        if (a>0) {
            real_areas.push_back(std::make_pair(i,a));
        }
    }
    real_areas.push_back(std::make_pair(coords_.size()-1,inf));
    
    size_t idx=0, curr=0; double min_area=0;
    
    size_t curr_len=real_areas.size();
    
    std::tie(idx,curr,min_area,curr_len) = find_and_drop_min(real_areas, curr_len);
    
    while ( min_area < inf) {
        //std::cout << idx << ", " << curr << ", " << min_area << ", " << curr_len << std::endl;
        //right
        size_t b = real_areas[idx-1].first;
        size_t c = real_areas[idx].first;
        size_t d = real_areas[idx+1].first;
        double right_area = triangle_area(coords_[b],coords_[c],coords_[d]);
        
        if (right_area < min_area) { right_area = min_area; }
        
        std::get<2>(coords_[c]) = right_area;
        real_areas[idx].second=right_area;
        
        if (idx>1) {
            //left
            size_t a = real_areas[idx-2].first;
            
            double left_area = triangle_area(coords_[a],coords_[b],coords_[c]);
        
            if (left_area < min_area) { left_area = min_area; }
            std::get<2>(coords_[b]) = left_area;
            real_areas[idx-1].second=left_area;
        }
        std::tie(idx,curr,min_area,curr_len) = find_and_drop_min(real_areas, curr_len);
    }
    
        
    std::sort(ordered_weights_.begin(),ordered_weights_.end(),[&](size_t l, size_t r) { return std::get<2>(coords_[l]) < std::get<2>(coords_[r]); });
    

}
        

std::vector<std::pair<double,double>> VWSimplify::by_threshold(double thres) {

    std::vector<std::pair<double,double>> result;
    for (const auto& c: coords_) {
        if (std::get<2>(c) >= thres) {
            result.push_back(std::make_pair(std::get<0>(c),std::get<1>(c)));
        }
    }
    
    return result;
}

std::shared_ptr<geos_geometry> simplify_linestring(std::shared_ptr<geos_geometry> geom, double thres) {
    if (!geom) { return geom; }
    
    auto cc=geom->get_coords();
    VWSimplify simp(cc);
    auto r = simp.by_threshold(thres);
    if (r.size()<2) {
        return nullptr;
    }
    return geos_geometry_create_linestring(geom->get_base(), r);
};

std::shared_ptr<geos_geometry> simplify_polygon(std::shared_ptr<geos_geometry> geom, double thres) {
    
    if (!geom) { return geom; }
    
    
    auto rings = geom->get_rings();
    
    if (rings.empty()) { return nullptr; }
    
    VWSimplify simp_outer(rings[0]);
    auto oo = simp_outer.by_threshold(thres);
    if (oo.size()<4) { return nullptr; }
    auto outer = geos_geometry_create_polygon(geom->get_base(), {oo});
    if (!outer->is_valid()) {
        outer = outer->buffer_zero();
    }
    
    if (rings.size()==1) {
        return outer;
    }
    
    std::vector<std::shared_ptr<geos_geometry>> interiors;
    for (size_t i=1; i < rings.size(); i++) {
        VWSimplify simp(rings[i]);
        auto r = simp.by_threshold(thres);
        
        if (r.size()>3) {
            auto g = geos_geometry_create_polygon(geom->get_base(), {r});
            if (!g->is_valid()) {
                g = g->buffer_zero();
            }
            interiors.push_back(g);
        }
        
    }
    if (interiors.empty()) {
        return outer;
    }
    auto interior = geos_geometry_create_collection(geom->get_base(), interiors)->unary_union();
    
    return outer->difference(interior);
}
    

std::shared_ptr<geos_geometry> simplify_geometry(std::shared_ptr<geos_geometry> in, double thres) {
    if (!in) { return in; }
    switch (in->type_id()) {
        case 0: return in;
        case 1: return in;
        case 2: return simplify_linestring(in, thres);
        case 3: return simplify_polygon(in, thres);
        case 4: return in;
    }
    std::vector<std::shared_ptr<geos_geometry>> result;
    for (auto g: in->get_geometries()) {
        if (g->type_id() == 2) {
            result.push_back(simplify_linestring(g, thres));
        } else if (g->type_id()==3) {
            auto r = simplify_polygon(g, thres);
            if (r) {
                result.push_back(r);
            }
        } else {
            auto r = simplify_geometry(g, thres);
            if (r) {
                result.push_back(r);
            }
        }
    }
    return geos_geometry_create_collection(in->get_base(), result);
}


        


namespace py = pybind11;
void export_vw(py::module& m) {
    py::class_<VWSimplify>(m, "VWSimplify")
        .def(py::init<std::vector<std::pair<double,double>>>())
        .def_property_readonly("coords", &VWSimplify::coords)
        .def_property_readonly( "ordered_weights",&VWSimplify::ordered_weights)
        .def("by_threshold", &VWSimplify::by_threshold)
    ;
    
    
    m.def("simplify_polygon_vw", simplify_polygon);
    m.def("simplify_linestring_vw", simplify_linestring);
    m.def("simplify_geometry_vw", simplify_geometry);
}
