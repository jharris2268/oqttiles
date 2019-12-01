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
#include <iterator>

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



std::shared_ptr<alt_tile_data> copy_alt_tile_data(std::shared_ptr<alt_tile_data> in) {
    
    auto result = std::make_shared<alt_tile_data>();
    result->tile_xyz = in->tile_xyz;
    for (auto o: in->objs) {
        
        result->objs.insert(o);
    }
    return result;
}



class AddBlocksTree {
    public:
        
    
        AddBlocksTree(std::shared_ptr<geos_geometry> geom_, alt_tile_data_vec_cb cb_, int64 maxdepth_) : geom(geom_), cb(cb_),maxdepth(maxdepth_) {}
        
        virtual ~AddBlocksTree() {}
        
        
                
        void call(std::shared_ptr<alt_tile_data> tl) {
            if (pending.size()>0) {
                if (tl && (tl->tile_xyz==curr_xyz)) {
                    //skip
                    
                } else {
                    auto tf = geometry_tiles_all_between(geom, curr_xyz, (tl ? tl->tile_xyz : xyz{-1,-1,-1}), maxdepth);
                    if (tf.empty()) {
                        //skip
                    }
                    std::cout << "\rcurr_xyz=" << curr_xyz << ", pending.size() = " << pending.size() << ", tf.size()=" << tf.size() << std::flush;
                    auto oo = std::make_shared<alt_tile_data_vec>();
                    oo->tile_xyz = curr_xyz;
                    oo->to_finish = tf;
                    for (auto& t: tf) {
                        auto it = pending.find(t);
                        if (it==pending.end()) {
                            //skip;
                        } else {
                            oo->tiles.push_back(it->second);
                            pending.erase(it);
                        }
                        
                    }
                    std::cout << "  pending.size() = " << pending.size() << ", oo->tiles.size()=" << oo->tiles.size() << std::flush;
                    cb(oo);
                }
            }
            if (!tl) {
                std::cout << "pending.size() = " << pending.size() << std::endl;
                cb(nullptr);
                return;
            }
            
            curr_xyz = tl->tile_xyz;
            for (auto o: tl->objs) {
                auto it = pending.find(o.first.first);
                if (it==pending.end()) {
                    
                    auto n = std::make_shared<alt_tile_data>();
                    n->tile_xyz = o.first.first;
                    n->objs.insert(o);
                    pending[o.first.first] = n;
                } else {
                    
                    auto jt = it->second->objs.find(o.first);
                    if (jt == it->second->objs.end()) {
                        it->second->objs.insert(o);
                    } else {
                        for (auto& x: o.second) {
                            jt->second.push_back(x);
                        }
                    }
                }                
                
            }
            
        }
        
    private:
        std::shared_ptr<geos_geometry> geom;
        alt_tile_data_vec_cb cb;
        int64 maxdepth;
        xyz curr_xyz;
        std::map<xyz,std::shared_ptr<alt_tile_data>> pending;
};
                
                
                


class AddBlocksTreeYY {
    public:
        
    
        AddBlocksTreeYY(std::shared_ptr<geos_geometry> geom_, alt_tile_data_vec_cb cb_, int64 maxdepth_) : geom(geom_), cb(cb_),maxdepth(maxdepth_) {}
        
        virtual ~AddBlocksTreeYY() {}
        
        
                
        void call(std::shared_ptr<alt_tile_data> tl) {
            if (curr) {
                if (tl && (tl->tile_xyz==curr->tile_xyz)) {
                    //skip
                    
                } else {
                    auto tf = geometry_tiles_all_between(geom, curr->tile_xyz, (tl ? tl->tile_xyz : xyz{-1,-1,-1}), maxdepth);
                    if (tf.empty()) {
                        //skip
                    }
                    
                    std::set<xyz> tfs(tf.begin(),tf.end());
                    
                    auto rr = std::make_shared<alt_tile_data>();
                    rr->tile_xyz = curr->tile_xyz;
                    for (auto t: curr->tiles) {
                        if (!t) {
                            
                            std::cout << "null t? tl xyz=" << (tl ? tl->tile_xyz : xyz{-1,-1,-1}) << " / curr xyz=" << curr->tile_xyz << std::endl;
                            throw std::domain_error("ZZZ");
                        } 
                        for (auto o: t->objs) {
                            if (tfs.count(o.first.first)>0) {
                                auto it = rr->objs.find(o.first);
                                if (it==rr->objs.end()) {
                                    rr->objs.insert(o);
                                } else {
                                    for (auto x: o.second) {
                                        it->second.push_back(x);
                                    }
                                }
                            }
                        }
                    }
                    auto oo = std::make_shared<alt_tile_data_vec>();
                    oo->tile_xyz = curr->tile_xyz;
                    oo->tiles.push_back(rr);
                    oo->to_finish = tf;
                    cb(oo);
                }
            }
                        
                            
                        
                    
            if (!tl) {
                cb(nullptr);
                return;
            }
            
            auto next = std::make_shared<alt_tile_data_vec>();
            next->tile_xyz = tl->tile_xyz;
            if (curr) {
                for (auto t: curr->tiles) {
                    if (!t) {
                        std::cout << "null t? tl xyz=" << tl->tile_xyz << " / curr xyz=" << curr->tile_xyz << std::endl;
                        throw std::domain_error("ZZ");
                    } else {
                        if (is_parent(t, tl)) {
                            next->tiles.push_back(t);
                        }
                    }
                }
            }
            next->tiles.push_back(tl);
            curr=next;
        }
            
            
            
    private:
        std::shared_ptr<geos_geometry> geom;
        alt_tile_data_vec_cb cb;
        int64 maxdepth;
        std::shared_ptr<alt_tile_data_vec> curr;
        
        
};



class AddBlocksTreeZZ {
    public:
        
    
        AddBlocksTreeZZ(std::shared_ptr<geos_geometry> geom_, alt_tile_data_vec_cb cb_, int64 maxdepth_) : geom(geom_), cb(cb_),maxdepth(maxdepth_) {
            
           
        }
        
        virtual ~AddBlocksTreeZZ() {}
        
        
                
        void call(std::shared_ptr<alt_tile_data> tl) {
            
            
            if (curr) {
                if (tl && (tl->tile_xyz==curr->tile_xyz)) {
                    //skip
                    
                } else {
                    auto tf = geometry_tiles_all_between(geom, curr->tile_xyz, (tl ? tl->tile_xyz : xyz{-1,-1,-1}), maxdepth);
                    if (tf.empty()) {
                        //skip
                    } else if (tf.size()<100) {
                        curr->to_finish = tf;
                        
                        cb(curr);
                        
                    } else {
                        for (size_t j=0; j < tf.size(); j+=50) {
                            size_t lj = std::min(tf.size(), j+50);
                            std::vector<xyz> tfc(tf.begin()+j, tf.begin()+lj);
                            auto curr_split = std::make_shared<alt_tile_data_vec>(curr->tile_xyz, curr->tiles, tfc);
                            cb(curr_split);
                            
                        }
                    }
                }
                
            }
            
            
            if (!tl) {
                cb(nullptr);
                return;
            }
            
            auto next = std::make_shared<alt_tile_data_vec>();
            next->tile_xyz = tl->tile_xyz;
            if (curr) {
                for (auto t: curr->tiles) {
                    if (!t) {
                        std::cout << "null t? tl xyz=" << tl->tile_xyz << " / curr xyz=" << curr->tile_xyz << std::endl;
                        throw std::domain_error("ZZ");
                    } else {
                        if (is_parent(t, tl)) {
                            next->tiles.push_back(t);
                        }
                    }
                }
            }
            next->tiles.push_back(tl);
            curr=next;
            malloc_trim(0);
        }
        
    private:
        std::shared_ptr<geos_geometry> geom;
        alt_tile_data_vec_cb cb;
        int64 maxdepth;
        std::shared_ptr<alt_tile_data_vec> curr;
        
        
};

            
std::function<void(std::shared_ptr<alt_tile_data>)> make_addblockstree_alt_cb(std::shared_ptr<geos_geometry> geom, alt_tile_data_vec_cb cb, int64 maxdepth) {
    
    auto abt=std::make_shared<AddBlocksTree>(geom,cb,maxdepth);
    return [abt](std::shared_ptr<alt_tile_data> tl) { abt->call(tl); };
}


class AddOtherFeats {
    public:
        AddOtherFeats(std::shared_ptr<alt_tile_data> otherfeatures_, alt_tile_data_vec_cb cb_)
            : otherfeatures(otherfeatures_), cb(cb_) {
            
            
            
            for (const auto& kv: otherfeatures->objs) {
                tabs.insert(kv.first.second);
                tiles[kv.first]=kv.second;
            }
            
        }
            
        
        void call(std::shared_ptr<alt_tile_data_vec> in) {
            if (!in) {
                if (!tiles.empty()) {
                    std::cout << "have " << tiles.size() << " otherfeats remaining" << std::endl;
                    
                    auto rr = std::make_shared<alt_tile_data_vec>();
                    rr->tile_xyz={0,0,0};
                    rr->to_finish.reserve(tiles.size());
                    for (const auto& kv: tiles) {
                        rr->to_finish.push_back(kv.first.first);
                    }
                    
                    
                    auto nn = std::make_shared<alt_tile_data>();
                    nn->tile_xyz = {0,0,0};nn->objs.swap(tiles);
                    
                    rr->tiles.push_back(nn);
                    
                    cb(rr);
                    
                }   
                
                cb(in);
                return;
            }
            
            auto nn = std::make_shared<alt_tile_data>();
            nn->tile_xyz = {0,0,0};
            for (const auto& tile: in->to_finish) {
                for (const auto& tab: tabs) {
                    auto key = std::make_pair(tile,tab);
                    //auto it = otherfeatures->objs.find(key);
                    auto it=tiles.find(key);
                    
                    if (it!=tiles.end()) {
                    
                        nn->objs[key] = it->second;
                        tiles.erase(it);
                    }
                    
                }
            }
            
            if (!nn->objs.empty()) {
                in->tiles.push_back(nn);
            }
            cb(in);
        }
            
    
    
    private:
        std::shared_ptr<alt_tile_data> otherfeatures;
        alt_tile_data_vec_cb cb;
        std::set<std::string> tabs;
        std::map<std::pair<xyz,std::string>,std::vector<alt_feature_data>> tiles;
};


            
alt_tile_data_vec_cb make_addotherfeatures(std::shared_ptr<alt_tile_data> otherfeatures, alt_tile_data_vec_cb cb) {
    if ((!otherfeatures) || (otherfeatures->objs.empty())) { return cb; }
    
    auto aof = std::make_shared<AddOtherFeats>(otherfeatures,cb);
    
    
    return [aof](std::shared_ptr<alt_tile_data_vec> in) {
        aof->call(in);
    };
}
        
        
    

    
void export_prepare_geometries(py::module& m) {
    
    py::class_<alt_tile_data_vec, std::shared_ptr<alt_tile_data_vec>>(m, "alt_tile_data_vec")
        .def(py::init<xyz,std::vector<std::shared_ptr<alt_tile_data>>,std::vector<xyz>>())
        .def_readonly("tile_xyz", &alt_tile_data_vec::tile_xyz)
        .def_readonly("tiles", &alt_tile_data_vec::tiles)
        .def_readonly("to_finish", &alt_tile_data_vec::to_finish)
    ;
    
    
    m.def("make_addblockstree_alt_cb", [](std::shared_ptr<geos_geometry> g, std::function<void(std::shared_ptr<alt_tile_data_vec>)> cb, int64 maxdepth) { return make_addblockstree_alt_cb(g,{cb},maxdepth); });
    m.def("make_addotherfeatures", &make_addotherfeatures);
    m.def("is_parent", [](std::shared_ptr<alt_tile_data> l, std::shared_ptr<alt_tile_data> r) { return is_parent(l,r); });
    
    //m.def("make_collectaltblocks_cb", [](std::shared_ptr<geos_geometry> g, std::function<void(std::shared_ptr<alt_tile_data_vec>)> cb, int64 maxdepth) { return make_collectaltblocks_cb(g,{cb},maxdepth); });
    
}


