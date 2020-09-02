
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>

#include <oqt/utils/threadedcallback.hpp>
#include <oqt/utils/multithreadedcallback.hpp>
#include <oqt/utils/splitcallback.hpp>

#include "oqttiles.hpp"
#include "geos_wrapper.hpp"
#include "maketiledata.hpp"
#include "prepare_geometries.hpp"
#include "mvt.hpp"

namespace py = pybind11;
using namespace std::literals::string_literals;

PYBIND11_DECLARE_HOLDER_TYPE(XX, std::shared_ptr<XX>);


class collect_packed_tiles {
    public:
        collect_packed_tiles(py::object cb_) : cb(cb_) {}
        
        virtual ~collect_packed_tiles() {}
        
        void cbwrap(std::shared_ptr<packed_tiles> data) {
            py::gil_scoped_acquire aq;
            py::dict pp;
            for (const auto& tl: *data) {
                pp[py::cast(tl.first)]=py::bytes(tl.second);
            }
            cb(pp);
        }
        
        void call(std::shared_ptr<packed_tiles> pt) {
            if (!pt) {
                if ((curr) && (curr->size()>0)) {
                    cbwrap(curr);
                }
                py::gil_scoped_acquire aq;
                cb(py::none());
                return;
            }
            
            if (!curr) {
                curr=pt;
            } else if (curr->size()>100) {
                
                cbwrap(curr);
                curr = pt;
            } else {
                for (const auto& t: *pt) {
                    curr->insert(t);
                }
            }
        }
    private:
        py::object cb;
        std::shared_ptr<packed_tiles> curr;
};


class WrapCb {
    public:
        WrapCb(const std::string& name_, std::function<void(oqt::PrimitiveBlockPtr)> cb_) : name(name_), cb(cb_), done(false) {}
        
        virtual ~WrapCb() {
            if (!done) {
                py::gil_scoped_release gr;
                std::cout << "~WrapCb() " << name << " not finished, call cb(nullptr)" << std::endl;
                cb(nullptr);
            }
        }
        
        std::string repr() { return name; }
        
        void call(std::vector<oqt::PrimitiveBlockPtr> bls) {
            
            if (done) {
                throw std::domain_error("already finished");
            }
            
            py::gil_scoped_release gr;
            
            if (bls.empty()) {
                done=true;
                cb(nullptr);
                return;
            }
            
            for (auto bl: bls) {
                if (!bl) {
                    throw std::domain_error("don't pass nulls");
                }
                cb(bl);
            }
        }
        
    private:
        std::string name;
        std::function<void(oqt::PrimitiveBlockPtr)> cb;
        bool done;
};
        
        
template <class T>
std::function<void(std::shared_ptr<T>)> wrap_callback(std::function<void(std::shared_ptr<T>)> cb, std::string what) {
    
    return [cb,what](std::shared_ptr<T> bl) {
        try {
            cb(bl);
        } catch(std::exception& ex) {
            std::cout << "?? " << what << ": " << ex.what() << std::endl;
        }
    };
}
        
    
        

std::shared_ptr<WrapCb> make_processall_alt_callback(
    std::shared_ptr<geos_base> gs, feature_spec features, std::map<std::string,extra_tags_spec> extra_tags,
    int64 max_tile, bool use_obj_tile, std::shared_ptr<geos_geometry> filter_poly, bool fix_geom, bool mergefeats,
    bool simplify_max, std::shared_ptr<alt_tile_data> otherfeatures, py::object cb) {

    size_t nt=4;
    if (!filter_poly) { throw std::domain_error("no filter_poly"); }
    
    auto cpt = std::make_shared<collect_packed_tiles>(cb);
    
    auto cpt2 = [cpt](std::shared_ptr<packed_tiles> pt) { cpt->call(pt); };
    auto cpt2w = wrap_callback<packed_tiles>(cpt2, "collect_packed_tiles");
    
    auto cpt_cb = oqt::multi_threaded_callback<packed_tiles>::make(cpt2w, nt);

    bbox filter_poly_bounds = filter_poly->bounds();
    double xb = (std::get<2>(filter_poly_bounds)-std::get<0>(filter_poly_bounds))*0.025;
    double yb = (std::get<3>(filter_poly_bounds)-std::get<1>(filter_poly_bounds))*0.025;
    bbox buffered_bounds = bbox{std::get<0>(filter_poly_bounds)-xb,std::get<1>(filter_poly_bounds)-yb,std::get<2>(filter_poly_bounds)+xb,std::get<3>(filter_poly_bounds)+yb};
        
    
    std::vector<alt_tile_data_vec_cb> mts;
    for (size_t i=0; i < nt; i++) {
        mts.push_back(oqt::threaded_callback<alt_tile_data_vec>::make(
            wrap_callback(
                make_merge_interim_tiles_callback(true, max_tile, fix_geom, mergefeats, cpt_cb[i]),
                "merge_interim_tiles "+std::to_string(i)
            )));
    }
    
    auto mt = oqt::split_callback<alt_tile_data_vec>::make(mts);
    if (otherfeatures) {
        
        mt = oqt::threaded_callback<alt_tile_data_vec>::make(wrap_callback(make_addotherfeatures(otherfeatures, mt),"addotherfeatures"));
    }
        
    
    auto abt = oqt::multi_threaded_callback<alt_tile_data>::make(wrap_callback(make_addblockstree_alt_cb(filter_poly, mt,max_tile),"addblockstree"), nt);
    
    std::vector<std::function<void(oqt::PrimitiveBlockPtr)>> pas;
    for(size_t i=0; i < nt; i++) {
        auto mtdc = make_maketiledata_alt_callback(features, extra_tags, max_tile, use_obj_tile, buffered_bounds,fix_geom,  simplify_max, abt[i]);
        
        auto mtdcw = wrap_callback(mtdc, "maketiledata_alt "+std::to_string(i));
        pas.push_back(oqt::threaded_callback<oqt::PrimitiveBlock>::make(mtdcw));
    }
        
   
    auto pa = oqt::split_callback<oqt::PrimitiveBlock>::make(pas);
    
    return std::make_shared<WrapCb>("processall_alt_callback", pa);
    
}



std::shared_ptr<WrapCb> make_processall_alt_callback_nt(
    std::shared_ptr<geos_base> gs, feature_spec features, std::map<std::string,extra_tags_spec> extra_tags,
    int64 max_tile, bool use_obj_tile, std::shared_ptr<geos_geometry> filter_poly, bool fix_geom, bool mergefeats, 
    bool simplify_max, std::shared_ptr<alt_tile_data> otherfeatures, py::object cb) {

    if (!filter_poly) { throw std::domain_error("no filter_poly"); }
    
    auto cpt = std::make_shared<collect_packed_tiles>(cb);
    
    auto cpt_cb = [cpt](std::shared_ptr<packed_tiles> pt) { cpt->call(pt); };

    bbox filter_poly_bounds = filter_poly->bounds();
    double xb = (std::get<2>(filter_poly_bounds)-std::get<0>(filter_poly_bounds))*0.025;
    double yb = (std::get<3>(filter_poly_bounds)-std::get<1>(filter_poly_bounds))*0.025;
    bbox buffered_bounds = bbox{std::get<0>(filter_poly_bounds)-xb,std::get<1>(filter_poly_bounds)-yb,std::get<2>(filter_poly_bounds)+xb,std::get<3>(filter_poly_bounds)+yb};
        
    
    alt_tile_data_vec_cb mt = make_merge_interim_tiles_callback(true, max_tile, fix_geom, mergefeats, cpt_cb);
    
    if (otherfeatures) {
        mt = make_addotherfeatures(otherfeatures, mt);
    }
    auto abt = make_addblockstree_alt_cb(filter_poly, mt,max_tile);
    
    auto pa = make_maketiledata_alt_callback(features, extra_tags, max_tile, use_obj_tile, buffered_bounds,fix_geom,  simplify_max, abt);
    
    return std::make_shared<WrapCb>("processall_alt_callback_nt", pa);
   
}

/*
std::shared_ptr<WrapCb> make_processall_groupalt_callback(
    std::shared_ptr<geos_base> gs, feature_spec features, std::map<std::string,extra_tags_spec> extra_tags,
    int64 max_tile, bool use_obj_tile, std::shared_ptr<geos_geometry> filter_poly, bool fix_geom, bool mergefeats, bool simplify_max, const std::vector<xyz>& qts,
    py::object cb) {

    size_t nt=4;
    if (!filter_poly) { throw std::domain_error("no filter_poly"); }
    
    auto cpt = std::make_shared<collect_packed_tiles>(cb);
    
    auto cpt_cb = oqt::threaded_callback<packed_tiles>::make([cpt](std::shared_ptr<packed_tiles> pt) { cpt->call(pt); }, nt);

    bbox filter_poly_bounds = filter_poly->bounds();
    double xb = (std::get<2>(filter_poly_bounds)-std::get<0>(filter_poly_bounds))*0.025;
    double yb = (std::get<3>(filter_poly_bounds)-std::get<1>(filter_poly_bounds))*0.025;
    bbox buffered_bounds = bbox{std::get<0>(filter_poly_bounds)-xb,std::get<1>(filter_poly_bounds)-yb,std::get<2>(filter_poly_bounds)+xb,std::get<3>(filter_poly_bounds)+yb};
        
    
    std::vector<alt_tile_data_vec_cb> mts;
    for (size_t i=0; i < nt; i++) {
        mts.push_back(oqt::threaded_callback<alt_tile_data_vec>::make(make_merge_interim_tiles_callback(true, max_tile, fix_geom, mergefeats, cpt_cb)));
    }
    
    
    auto abt = oqt::multi_threaded_callback<alt_tile_data_vec>::make(make_collectaltblocks_cb(filter_poly,mts,max_tile), nt);
    
    std::vector<std::function<void(oqt::PrimitiveBlockPtr)>> pas;
    for(size_t i=0; i < nt; i++) {
        auto mtdc = make_maketiledata_groupalt_callback(features, extra_tags, max_tile, use_obj_tile, buffered_bounds,fix_geom,  simplify_max, qts, abt[i]);
        pas.push_back(oqt::threaded_callback<oqt::PrimitiveBlock>::make(mtdc));
    }
        
    
    auto pa = oqt::split_callback<oqt::PrimitiveBlock>::make(pas);
    
    return std::make_shared<WrapCb>("processall_alt_callback", pa);
   
}
*/

    
    
    


void export_oqttile(py::module& m) {
    

    m.def("tile_bound", &tile_bound);
    
    

    m.def("bbox_overlaps", &bbox_overlaps);
    
    m.def("make_processall_alt_callback", &make_processall_alt_callback);
    m.def("make_processall_alt_callback_nt", &make_processall_alt_callback_nt);
    //m.def("make_processall_groupalt_callback", &make_processall_groupalt_callback);
    
    py::class_<WrapCb,std::shared_ptr<WrapCb>>(m,"WrapCb")
        .def("__repr__", &WrapCb::repr)
        .def("__call__", &WrapCb::call)
    ;
}



void export_geos_geometry(py::module& m);
void export_prepare_geometries(py::module& m);
void export_mvt(py::module& m);
void export_maketiledata(py::module& m);

void export_vw(py::module& m);

PYBIND11_MODULE(_oqttiles, m) {
    
    
    //py::module m("_oqttiles", "oqttiles");
    export_oqttile(m);
    
    export_prepare_geometries(m);
    export_geos_geometry(m);
    export_mvt(m);
    export_maketiledata(m);
    export_vw(m);
    //return m.ptr();
}



