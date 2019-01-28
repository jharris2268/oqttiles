#include "prepare_geometries.hpp"


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <oqt/elements/geometry.hpp>
#include <oqt/elements/block.hpp>
#include <oqt/elements/quadtree.hpp>
#include <oqt/geometry/elements/point.hpp>
#include <oqt/geometry/elements/linestring.hpp>
#include <oqt/geometry/elements/simplepolygon.hpp>
#include <oqt/geometry/elements/complicatedpolygon.hpp>
#include <oqt/utils/timing.hpp>
#include <oqt/utils/operatingsystem.hpp>
#include <malloc.h>
#include <fstream>

#include "maketiledata.hpp"
using namespace std::literals::string_literals;
PYBIND11_DECLARE_HOLDER_TYPE(XX, std::shared_ptr<XX>);

namespace py = pybind11;

template <class td_type> int64 tx(std::shared_ptr<td_type> l) { return std::get<0>(l->tile_xyz); }
template <class td_type> int64 ty(std::shared_ptr<td_type> l) { return std::get<1>(l->tile_xyz); }
template <class td_type> int64 tz(std::shared_ptr<td_type> l) { return std::get<2>(l->tile_xyz); }

template <class td_type>
bool is_parent_int(std::shared_ptr<td_type> l, std::shared_ptr<td_type> r) {
    if ((!l) || (!r)) { throw std::domain_error("null tile_data!"); }
        
    if (tz(l) > tz(r)) { return false; }
    
    if (tz(l) == tz(r)) {
        return (tx(l) == tx(r)) && ( ty(l) == ty(r));
    }
    
    int64 p = 1ll << (tz(r) - tz(l));
    int64 rx = tx(r) / p;
    int64 ry = ty(r) / p;
    
    return (tx(l)==rx) && (ty(l)==ry);
}

//bool is_parent(std::shared_ptr<tile_data> l, std::shared_ptr<tile_data> r) { return is_parent_int(l,r); }
bool is_parent(std::shared_ptr<alt_tile_data> l, std::shared_ptr<alt_tile_data> r) { return is_parent_int(l,r); }


template <class td_type>
class AddBlocksTree {
    public:
        typedef typename std::function<void(std::shared_ptr<data_vec<td_type>>)> callback_type;
    
    
        AddBlocksTree(std::shared_ptr<geos_geometry> geom_, std::vector<callback_type> cb_, int64 maxdepth_) : geom(geom_), cb(cb_),maxdepth(maxdepth_) {
            //curr = std::make_shared<tile_data_vec>();
            i=0;
        }
        virtual ~AddBlocksTree() {}
        void call(std::shared_ptr<td_type> tl) {
            
            
            if (curr) {
                if (tl && (tl->tile_xyz==curr->tile_xyz)) {
                    //skip
                    //std::cout << "duplicate tile " << tup_str(curr->tile_xyz) << std::endl;
                } else {
                    auto tf = geometry_tiles_all_between(geom, curr->tile_xyz, (tl ? tl->tile_xyz : xyz{-1,-1,-1}), maxdepth);
                    if (tf.empty()) {
                        //std::cout << "skip tile_data_vec " << tup_str(curr->tile_xyz) << " "<< curr->tiles.size() << " tiles, none to finish" << std::endl;
                        //pass
                    } else if (tf.size()<100) {
                        curr->to_finish = tf;
                        
                        cb[i%cb.size()](curr);
                        i++;
                    } else {
                        //std::cout << "split tile_data_vec " << tup_str(curr->tile_xyz) << " " << curr->tiles.size() << " tiles, " << tf.size() << " to finish" << std::endl;
                        for (size_t j=0; j < tf.size(); j+=50) {
                            size_t lj = std::min(tf.size(), j+50);
                            std::vector<xyz> tfc(tf.begin()+j, tf.begin()+lj);
                            auto curr_split = std::make_shared<data_vec<td_type>>(curr->tile_xyz, curr->tiles, tfc);
                            cb[i%cb.size()](curr_split);
                            i++;
                        }
                    }
                }
                
            }
            
            
            if (!tl) {
                for (auto c:cb) {
                    c(nullptr);
                }
                return;
            }
            
            auto next = std::make_shared<data_vec<td_type>>();
            next->tile_xyz = tl->tile_xyz;
            if (curr) {
                for (auto t: curr->tiles) {
                    if (is_parent(t, tl)) {
                        next->tiles.push_back(t);
                    }
                }
            }
            next->tiles.push_back(tl);
            curr=next;
            malloc_trim(0);
        }
        
    private:
        std::shared_ptr<geos_geometry> geom;
        std::vector<callback_type> cb;
        int64 maxdepth;
        std::shared_ptr<data_vec<td_type>> curr;
        size_t i=0;
};

            
std::function<void(std::shared_ptr<alt_tile_data>)> make_addblockstree_alt_cb(std::shared_ptr<geos_geometry> geom, std::vector<alt_tile_data_vec_cb> cb, int64 maxdepth) {
    
    auto abt=std::make_shared<AddBlocksTree<alt_tile_data>>(geom,cb,maxdepth);
    return [abt](std::shared_ptr<alt_tile_data> tl) { abt->call(tl); };
}



class CollectAltBlocks {
    public:
        CollectAltBlocks(std::shared_ptr<geos_geometry> geom_, std::vector<alt_tile_data_vec_cb> cb_, int64 max_depth_ )
            : geom(geom_), cb(cb_), max_depth(max_depth_), i(0), outss("collectaltblocks.log", std::ios::out) { pid=getpid(); }
        
        
        virtual ~CollectAltBlocks() {}
        
        void call(std::shared_ptr<alt_tile_data_vec> vals) {
            if (!vals) {
                finish(xyz{-1,-1,-1});
                
                
                outss << ts << " CollectAltBlocks: finished: " << temps.size() << " blocks left [" << oqt::getmem(pid) << "]" << std::endl;
                for (auto c:cb) {
                    c(nullptr);
                }
                return;
            }
            
            finish(vals->tile_xyz);
            
            for (const auto& xx: vals->tiles) {
                int64 q=oqt::quadtree::from_tuple(std::get<0>(xx->tile_xyz),std::get<1>(xx->tile_xyz),std::get<2>(xx->tile_xyz));
                temps.insert(std::make_pair(q, xx));
            }
            
            outss << ts << " CollectAltBlocks: " << temps.size() << " blocks";
            if (!temps.empty()) {
                outss << " " << tup_str(oqt::quadtree::tuple(temps.begin()->first)) << " to " << tup_str(oqt::quadtree::tuple(temps.rbegin()->first));
            }
            outss <<  " [" << oqt::getmem(pid) << "]" << std::endl;
        }
        
        void prep_finish(std::shared_ptr<alt_tile_data_vec> result, xyz last) {
            int64 lq = (std::get<2>(last)==-1) ? -1 : oqt::quadtree::from_tuple(std::get<0>(last),std::get<1>(last),std::get<2>(last));
            while (!temps.empty()) {
                auto it = temps.begin();
                if (!( (lq==-1) || (it->first<lq))) {
                    return;
                }
                result->tiles.push_back(it->second);
                temps.erase(it);
            }
        }
        
        void finish(xyz last) {
            if (temps.empty()) {
                return;
            }
            
            auto result = std::make_shared<alt_tile_data_vec>();
            
            result->tile_xyz = oqt::quadtree::tuple(temps.begin()->first);
            result->to_finish = geometry_tiles_all_between(geom, result->tile_xyz, last, max_depth);
            if (result->to_finish.empty()) { return; }
            prep_finish(result, last);
            malloc_trim(0);
            
            if (result->tiles.empty()) { return; }
            outss << ts << " CollectAltBlocks: " << tup_str(result->tile_xyz) << " " << result->to_finish.size() << " to finish, " << result->tiles.size() << " blocks, " << temps.size() << " remaining [" << oqt::getmem(pid) << "]" << std::endl;
            
            
            cb[i%cb.size()](result);
            i++;
        }
    private:
        std::shared_ptr<geos_geometry> geom;
        std::vector<alt_tile_data_vec_cb> cb;
        int64 max_depth;
        size_t i;
        std::ofstream outss;
        oqt::TimeSingle ts;
        std::multimap<int64, std::shared_ptr<alt_tile_data>> temps;
        size_t pid;
};
            
std::function<void(std::shared_ptr<alt_tile_data_vec>)> make_collectaltblocks_cb(std::shared_ptr<geos_geometry> geom, std::vector<alt_tile_data_vec_cb> cb, int64 maxdepth) {
    
    auto abt=std::make_shared<CollectAltBlocks>(geom,cb,maxdepth);
    return [abt](std::shared_ptr<alt_tile_data_vec> tl) { abt->call(tl); };
}            
            
            
            


void export_prepare_geometries(py::module& m) {
    
    py::class_<alt_tile_data_vec, std::shared_ptr<alt_tile_data_vec>>(m, "alt_tile_data_vec")
        .def(py::init<xyz,std::vector<std::shared_ptr<alt_tile_data>>,std::vector<xyz>>())
        .def_readonly("tile_xyz", &alt_tile_data_vec::tile_xyz)
        .def_readonly("tiles", &alt_tile_data_vec::tiles)
        .def_readonly("to_finish", &alt_tile_data_vec::to_finish)
    ;
    
    
    m.def("make_addblockstree_alt_cb", [](std::shared_ptr<geos_geometry> g, std::function<void(std::shared_ptr<alt_tile_data_vec>)> cb, int64 maxdepth) { return make_addblockstree_alt_cb(g,{cb},maxdepth); });
    
    m.def("is_parent", [](std::shared_ptr<alt_tile_data> l, std::shared_ptr<alt_tile_data> r) { return is_parent(l,r); });
    
    m.def("make_collectaltblocks_cb", [](std::shared_ptr<geos_geometry> g, std::function<void(std::shared_ptr<alt_tile_data_vec>)> cb, int64 maxdepth) { return make_collectaltblocks_cb(g,{cb},maxdepth); });
    
}


