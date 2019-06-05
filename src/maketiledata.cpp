#include "maketiledata.hpp"


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
#include <oqt/utils/logger.hpp>

#include "mvt.hpp"



using namespace std::literals::string_literals;
PYBIND11_DECLARE_HOLDER_TYPE(XX, std::shared_ptr<XX>);

namespace py = pybind11;



int geom_type(const std::shared_ptr<oqt::BaseGeometry>& geom) {
    if (geom->Type()==oqt::ElementType::Point) { return 0; }
    if (geom->Type()==oqt::ElementType::Linestring) { return 1; }
    if (geom->Type()==oqt::ElementType::SimplePolygon) { return 2; }
    if (geom->Type()==oqt::ElementType::ComplicatedPolygon) { return 2; }
    throw std::domain_error("wrong elementtype??");
    return -1;
}

bool has_tag(const std::vector<oqt::Tag>& tv, const oqt::Tag& tg) {
    for (const auto& t: tv) {
        if (t.key==tg.key) {
            return true;
        }
    }
    return false;
}
std::string tup_str(const xyz& t) {
    std::stringstream ss;
    ss << "[" << std::setw(6) << std::get<0>(t) << ", " << std::setw(6) << std::get<1>(t) <<  ", " << std::setw(2) << std::get<2>(t) << "]";
    return ss.str();
}

bool same_geomtype(std::shared_ptr<geos_geometry> left, std::shared_ptr<geos_geometry> right) {
    int li = left->type_id();
    int ri = right->type_id();
    if (li==ri) { return true; }
    
    if ((li==4) && (ri==0)) { return true; }
    if ((li==5) && (ri==1)) { return true; }
    if ((li==6) && (ri==3)) { return true; }
    if ((li==0) && (ri==4)) { return true; }
    if ((li==1) && (ri==5)) { return true; }
    if ((li==3) && (ri==6)) { return true; }
    return false;
}

std::shared_ptr<geos_geometry> as_type(std::shared_ptr<geos_geometry> geom, int64 ty) {
    if (ty<0) { return geom; }
    
    if (ty==0) {
        if ((geom->type_id() == 0) || (geom->type_id()==4)) { return geom; }
        
        if ((geom->type_id() == 1) || (geom->type_id()==3)) { return geom->centroid(); }
        if ((geom->type_id() == 5) || (geom->type_id()==6)) {
            std::vector<std::shared_ptr<geos_geometry>> gg;
            for (auto g: geom->get_geometries()) {
                gg.push_back(g->centroid());
            }
            return geos_geometry_create_collection(geom->get_base(), gg);
        } else {
            return nullptr;
        }
    }
    
    if (ty==1) {
        if ((geom->type_id() == 0) || (geom->type_id()==4)) { return nullptr; }
        if ((geom->type_id() == 1) || (geom->type_id()==5)) { return geom; }
        if ((geom->type_id() == 3) || (geom->type_id()==6)) { return geom->boundary(); }
        return nullptr;
    }
    
    if (ty==2) {
        if ((geom->type_id() == 0) || (geom->type_id()==4)) { return nullptr; }
        if ((geom->type_id() == 1) || (geom->type_id()==5)) { return nullptr; }
        if ((geom->type_id() == 3) || (geom->type_id()==6)) { return geom; }
        return nullptr;
    }
    throw std::domain_error("??");
    return nullptr;
}

bool bbox_contains(const bbox& outer, const bbox& inner) {
    if (std::get<0>(outer)>std::get<0>(inner)) { return false; }
    if (std::get<1>(outer)>std::get<1>(inner)) { return false; }
    if (std::get<2>(outer)<std::get<2>(inner)) { return false; }
    if (std::get<3>(outer)<std::get<3>(inner)) { return false; }
    return true;
}
std::shared_ptr<geos_geometry> transform_tile(std::shared_ptr<geos_geometry> geom, const bbox& bx) {
    
    double x0 = -1*std::get<0>(bx);
    double y0 = -1*std::get<1>(bx);
    double sc = 256./(std::get<2>(bx)-std::get<0>(bx));
    return  geos_geometry_transform_translate_scale(geom, x0, y0, sc);
}


std::shared_ptr<geos_geometry> clip_geom_tile(std::shared_ptr<geos_geometry> gsimp, xyz tt, bool checkfull, bool transform, bool fix_clipped) {
    auto bx0 = tile_bound(std::get<0>(tt),std::get<1>(tt),std::get<2>(tt),0.);
    
    auto obj_bx = gsimp->bounds();
    if (bbox_contains(bx0, obj_bx)) {
        
        if (transform) {
            return transform_tile(gsimp, bx0);
        }
        return gsimp;
    }
    if (gsimp->type_id()==0) { return nullptr; }
    if (gsimp->type_id()==4) {
        auto gg=gsimp->get_geometries();
        std::vector<std::shared_ptr<geos_geometry>> rr;
        for (auto g: gg){
            auto r=clip_geom_tile(g,tt,checkfull, transform, fix_clipped);
            if (r) { rr.push_back(r); }
        }
        if (rr.empty()) { return nullptr; }
        if (rr.size()==1) { return rr[0]; }
        return geos_geometry_create_collection(gsimp->get_base(), rr);
    }
    
    
    auto gc0 = gsimp->clip(std::get<0>(bx0),std::get<1>(bx0),std::get<2>(bx0),std::get<3>(bx0));
    if (!gc0) { return nullptr; }
    if (gc0->is_empty()) { return nullptr; }
    
    if (!checkfull) {
        if (transform) {
            return transform_tile(gc0, bx0);
        }
        return gc0;
    }
    
    auto bx1 = tile_bound(std::get<0>(tt),std::get<1>(tt),std::get<2>(tt),0.05);
    auto gc1 = gsimp->clip(std::get<0>(bx1),std::get<1>(bx1),std::get<2>(bx1),std::get<3>(bx1));
    if (!gc1 || gc1->is_empty()) {
        gc1=gc0;
    }
    if (gc1->type_id()==7) {
        //std::cout << "filter_collection " << " clipped tile " << tup_str(tt) << std::endl;
        gc1=filter_collection(gsimp->type_id(), gc1);
        if (gc1->type_id()==7) {
            //std::cout << " failed?? gsimp->type_()=" << gc1->type_() << std::endl;
            return nullptr;
        }
        
    }
    if (!same_geomtype(gsimp, gc1)) {
        std::cout << "skip geometry [ " << gsimp->type_() << "!=" << gc1->type_() << "] clipped tile " << tup_str(tt) << std::endl;
        return nullptr;
    }
    
    if (transform) {
        auto gt=transform_tile(gc1, bx0);
        if (fix_clipped && (!gt->is_valid())) {
            gt = gt->buffer_zero();
        }
        return gt;
    }
    if (fix_clipped && (!gc1->is_valid())) {
        gc1=gc1->buffer_zero();
    }
    return gc1;
}    

std::shared_ptr<geos_geometry> simplify_geom_check(std::shared_ptr<geos_geometry>geom, int64 zoom) {
    double simp = zoom_scale(zoom)/4.0;
    
    auto gsimp = geom->simplify(simp);

    if (!gsimp) {
        return nullptr;
    }
    if (gsimp->type_id() == 7) {
        if (gsimp->is_empty()) {
            return nullptr;
        }
        std::cout << "filter_collection " << /*tup_str(basetile) << " " << geom->Type() << " " << geom->Id() <<*/ " simp zoom " << zoom << std::endl;
        gsimp = filter_collection(geom->type_id(), gsimp);
        if (gsimp->type_id()==7) {
            std::cout << " failed?? gsimp->type_()=" << gsimp->type_() << std::endl;
        }
        
        if (gsimp->is_empty()) {
            return nullptr;
        }
    }
    
    if (!same_geomtype(geom, gsimp)) {
        std::cout << "skip geometry [ " << geom->type_() << "!=" << gsimp->type_() << "] " << /*tup_str(basetile) << " " << geom->Type() << " " << geom->Id() <<*/ " simp zoom " << zoom << std::endl;
        return nullptr;
    }
    if (!gsimp->is_valid()) {
        gsimp = gsimp->buffer_zero();
    }
    /*if (!gsimp->is_valid()) {
        std::cout << "simplify still invalid? try buffer(0.0001):" << std::flush;
        gsimp = gsimp->buffer(0.0001,16,1,1,0);
        std::cout << " is_valid" << (gsimp->is_valid() ? "t"s : "f"s) << std::endl;
    }*/
    
    
    gsimp->bounds();
    
    return gsimp;
}

std::shared_ptr<geos_geometry> make_geometry(std::shared_ptr<geos_base> gst, const bbox& filter_box, std::shared_ptr<oqt::BaseGeometry> geom, int64 geom_type) {
    auto orig_geom = geos_geometry_create_wkb(gst, geom->Wkb(true, false));
    if (!orig_geom) { return nullptr; }
    
    if (!orig_geom->is_valid()) {
        orig_geom=orig_geom->buffer_zero();
    }
    if (!orig_geom->is_valid()) {
        return nullptr;
    }

    orig_geom = as_type(orig_geom, geom_type);
    if (!orig_geom) { return nullptr; }
    
    orig_geom = orig_geom->clip(std::get<0>(filter_box), std::get<1>(filter_box), std::get<2>(filter_box), std::get<3>(filter_box));
    if (!orig_geom) { return nullptr; }
    if (orig_geom->is_empty()) { return nullptr; }
    return orig_geom;
}

std::map<xyz,std::shared_ptr<geos_geometry>> prep_geometries_noclip(std::shared_ptr<geos_geometry> geom, int64 minzoom, int64 max_tile, bool simp_max) {
    std::map<xyz, std::shared_ptr<geos_geometry>> geom_simp;
    int64 mz=(minzoom<max_tile) ? minzoom : max_tile;
    
    if (simp_max) {
        geom_simp[xyz{-1,-1,max_tile}] = simplify_geom_check(geom, max_tile);
    } else {
        geom_simp[xyz{-1,-1,max_tile}] = geom;
    }
        
    for (int64 zoom = mz; zoom<max_tile; zoom++) {
        auto gsimp = simplify_geom_check(geom, zoom);
        if (gsimp) {
            geom_simp[xyz{-1,-1,zoom}] = gsimp;
        }
    }
    
    return geom_simp;
}

bool is_box(std::shared_ptr<geos_geometry> geom) {
    if (geom->type_id()!=3) { return false; }
    auto rr = geom->get_rings();
    if (rr.size()!=1) { return false; }
    if (rr[0].size()!=5) { return false; }
    
    return geom->contains(geos_geometry_create_box(geom->get_base(), 0, 0, 256, 256));
}

size_t add_all_children(xyz basetile, std::shared_ptr<geos_geometry> geom, int64 min_zoom, int64 max_tile, const std::function<void(xyz t, std::shared_ptr<geos_geometry> g)>& cb) {
    xyz last= next_branch(basetile);
    auto tt = next_tile(basetile, max_tile);
    //std::cout << "call add_all_children " << tup_str(basetile) << " max=" << max_tile << " " << tup_str(tt) << " to " << tup_str(last) << std::endl;
    if (tt==last) { return 0; }
    
    size_t r=0;
    for ( ; tile_lessthan(tt, last); tt = next_tile(tt,max_tile)) {
        if (std::get<2>(tt)>=min_zoom) {
            cb(tt,geom);
        }
        r++;
    }
    return r;
}
    

void prep_geometries_clip_cb(std::shared_ptr<geos_geometry> geom, xyz basetile, int64 minzoom, int64 max_tile, bool fix_geoms, bool check_box, const std::function<void(xyz t, std::shared_ptr<geos_geometry> g)>& cb, bool simplify_max) {
    
    auto geom_simp_temp = prep_geometries_noclip(geom,  minzoom, max_tile, simplify_max);
    
    if (std::get<2>(basetile) > minzoom) {
        for (int64 zoom=minzoom; zoom < std::get<2>(basetile); zoom++) {
            
            xyz tt = tile_round(basetile, zoom);
            if (geom_simp_temp.count(xyz{-1,-1,std::get<2>(tt)})>0) {
            
                std::shared_ptr<geos_geometry> gg = geom_simp_temp.at(xyz{-1,-1,std::get<2>(tt)});
                
                if (!gg) {
                    throw std::domain_error("???? not geom for "+tup_str(tt));
                }
                auto gc = clip_geom_tile(gg, tt, true, true, fix_geoms);
                if (gc) {
                    
                    cb(tt,gc);
                }
            }
        }
    }
        
    xyz tt = basetile;
    xyz lt = next_branch(tt);
    
    try {
        while (tile_lessthan(tt,lt)) {
            std::shared_ptr<geos_geometry> gg;
            xyz k{-1,-1,std::get<2>(tt)};
            
            if (std::get<2>(tt) < minzoom) {
                gg = geom_simp_temp.begin()->second;
            } else if (geom_simp_temp.count(k)>0) {
                gg = geom_simp_temp.at(k);
            }
            if (!gg) {
                throw std::domain_error("???? not geom for "+tup_str(tt));
            }

            auto gc = clip_geom_tile(gg, tt,true, true, fix_geoms);
            
            auto x=tt;
            if (gc) {
                if (geom_simp_temp.count(k)>0) {
                    cb(tt,gc);
                }
                if (check_box && is_box(gc) && (std::get<2>(tt)<max_tile)) {
                    /*size_t n =*/ add_all_children(tt, gc, minzoom, max_tile, cb);
                    //std::cout << tup_str(tt) << " added " << n << " children" << std::endl;
                    tt = next_branch(tt);
                } else {
                        
                    tt=next_tile(tt,max_tile);
                }
            } else {
                tt=next_branch(tt);
            }
            if (std::get<0>(tt) < 0) {
                if (tt!=lt) { 
                    std::cout << "?? " << x << " next tile = " << tt << std::endl;
                }
            }
            
        }
    } catch(std::exception& ex) {
        std::cout << "failed at tile " << tt << std::endl;
        throw ex;
    }

}
std::map<xyz,std::shared_ptr<geos_geometry>> prep_geometries_clip(std::shared_ptr<geos_geometry> geom, xyz basetile, int64 minzoom, int64 max_tile, bool fix_geoms, bool check_box, bool simplify_max) {
        
        
    std::map<xyz, std::shared_ptr<geos_geometry>> geom_simp;        
    prep_geometries_clip_cb(geom, basetile, minzoom, max_tile, fix_geoms, check_box, [&geom_simp](xyz t, std::shared_ptr<geos_geometry> g) { geom_simp[t]=g; },simplify_max);
    return geom_simp;
}



class FeatureProperties {

    public:
        FeatureProperties(
            feature_spec features_,
            std::map<std::string,extra_tags_spec> extra_tags_,
            int64 max_tile__)
            
            : features(features_), extra_tags(extra_tags_), max_tile_(max_tile__) {
                
                for (const auto& f: features) {
                    for (const auto& t: f.second.second) {
                        
                        if (extra_tags.count(t)==0) {
                            throw std::domain_error("unrecognised table "+t);
                        }
                    }
                }
            }
            
        virtual ~FeatureProperties() {}
        
        int64 max_tile() const { return max_tile_; }
        
        void check_feature(std::map<std::string, std::vector<oqt::Tag>>& tabs, int ty, const oqt::Tag& tg, int64 zoom) const {
            
            auto it = features.find(std::make_tuple(ty, tg.key, tg.val));
            if (it!=features.end()) {
                if ((zoom==-1) || (zoom >= it->second.first)) {
                    
                    for (const auto& tab: it->second.second) {
                        tabs[tab].push_back(tg);
                    }
                }
            }
            
            it = features.find(std::make_tuple(ty, tg.key, "*"s));
            if (it!=features.end()) {
                if ((zoom==-1) || (zoom >= it->second.first)) {
                    for (const auto& tab: it->second.second) {
                        tabs[tab].push_back(tg);
                    }
                }
            }
            
            
        }
        bool check_feat_tags(const std::vector<oqt::Tag>& src, const std::pair<std::string,std::string>& test) const {
            for (const auto& t: src) {
                if ((t.key==test.first) && ((test.second=="*"s) || (t.val==test.second))) { return true; }
            }
            return false;
        }
        bool check_feat_tags_all(const std::vector<oqt::Tag>& src, const std::vector<std::pair<std::string,std::string>>& tests) const {
            for (const auto& test: tests) {
                if (!check_feat_tags(src, test)) { return false;}
            }
            return true;
        }
        void check_extra_tags(const std::string& tb, std::vector<oqt::Tag>& tgs, const oqt::Tag& tg, int64 zoom) const {
            if (has_tag(tgs, tg)) { return; }
            if (zoom==-1) {
                tgs.push_back(tg);
                return;
            }
            
            auto extra_spec = extra_tags.find(tb);
            if (extra_spec==extra_tags.end()) { return; }
            
            
            
            for (const auto& sp: extra_spec->second.second) {
                
                if (std::get<1>(sp) > zoom) { continue; }
                
                if ((std::get<2>(sp) != "*") && (std::get<2>(sp) != tg.key)) { continue; }
                
                if (!std::get<0>(sp).empty()) {
                    if (!check_feat_tags_all(tgs, std::get<0>(sp))) {
                        continue;
                    }
                }
                
                tgs.push_back(tg);
            }
        }
                
        feature_tables prep_feature_table(std::shared_ptr<oqt::BaseGeometry> geom, int ty) const {
            feature_tables tabs;
            auto mm = prep_feature_table_entries(geom,max_tile()==14 ? -1 : max_tile(),ty);
            if (mm.empty()) { return {}; }
            
            tabs[max_tile()]=mm;
            int64 minzoom=geom->MinZoom();
            
            if (minzoom < max_tile()) {
                for (int64 zoom = minzoom; zoom < max_tile(); zoom++) {
                    auto e = prep_feature_table_entries(geom, zoom, ty);
                    if (!e.empty()) {
                        tabs[zoom] = e;
                    }
                }
            }
            
            return tabs;
        }
        
        int table_geom_type(const std::string& t) const {
            if (extra_tags.count(t)==0) { throw std::domain_error("unknown table "+t); }
            return extra_tags.at(t).first;
        }
        std::vector<feature_table_entry> prep_feature_table_entries(std::shared_ptr<oqt::BaseGeometry> geom, int64 zoom, int ty) const {
            
            std::map<std::string,std::vector<oqt::Tag>> tabs;
            int gt=geom_type(geom);
            for (const auto& tg: geom->Tags()) {
                check_feature(tabs, gt, tg, zoom);
                    
            }
            
            std::vector<feature_table_entry> result;
            for (auto& t: tabs) {
                
                //if ( (ty>=0) && (t.second.first!=ty)) { continue; }
                if ((ty>=0) && (ty!=table_geom_type(t.first))) { continue; }
                
                
                
                
                for (const auto& tg: geom->Tags()) {
                    check_extra_tags(t.first, t.second, tg, zoom);
                }
                   
                feature_table_entry entry;
                std::get<0>(entry) = t.first;
                
                std::stringstream key_strm;
                
                std::sort(t.second.begin(),t.second.end(),[](const oqt::Tag&l, const oqt::Tag& r) { return l.key<r.key; });
                
                
                property_map& pm = std::get<2>(entry);
                
                for (const auto& tg: t.second) {
                    picojson::value v(tg.val);
                    if (zoom>=0) {
                        key_strm << tg.key << "=" << v.to_str() << ";";
                    }
                    pm[tg.key]=v;
                }
                if (geom->Type() == oqt::ElementType::Point) {
                    auto pt = std::dynamic_pointer_cast<oqt::geometry::Point>(geom);
                    if (!pt) { throw std::domain_error("not a Point"); }
                
                    if (pt->Layer() != 0) {
                        if (zoom>=0) { key_strm << "layer=" << pt->Layer(); }
                        pm["layer"] = picojson::value(pt->Layer());
                    }
                }
                
                if (geom->Type() == oqt::ElementType::Linestring) {
                    auto ls = std::dynamic_pointer_cast<oqt::geometry::Linestring>(geom);
                    if (!ls) { throw std::domain_error("not a Linestring"); }
                    
                    
                    if (ls->ZOrder() != 0) {
                        if (zoom>=0) { key_strm << "z_order=" << ls->ZOrder(); }
                        pm["z_order"] = picojson::value(ls->ZOrder());
                    }
                    
                    if (ls->Layer() != 0) {
                        if (zoom>=0) { key_strm << "layer=" << ls->Layer(); }
                        pm["layer"] = picojson::value(ls->Layer());
                    }
                    pm["way_length"] = picojson::value(ls->Length());
                }
                
                if (geom->Type() == oqt::ElementType::SimplePolygon) {
                    auto sp = std::dynamic_pointer_cast<oqt::geometry::SimplePolygon>(geom);
                    if (!sp) { throw std::domain_error("not a SimplePolygon"); }
                    
                    
                    if (sp->ZOrder() != 0) {
                        if (zoom>=0) { key_strm << "z_order=" << sp->Layer(); }
                        pm["z_order"] = picojson::value(sp->ZOrder());
                    }
                    
                    if (sp->Layer() != 0) {
                        if (zoom>=0) { key_strm << "layer=" << sp->Layer(); }
                        pm["layer"] = picojson::value(sp->Layer());
                    }
                    pm["way_area"] = picojson::value(sp->Area());
                }
                
                if (geom->Type() == oqt::ElementType::ComplicatedPolygon) {
                    auto sp = std::dynamic_pointer_cast<oqt::geometry::ComplicatedPolygon>(geom);
                    if (!sp) { throw std::domain_error("not a ComplicatedPolygon"); }
                    
                    
                    if (sp->ZOrder() != 0) {
                        if (zoom>=0) { key_strm << "z_order=" << sp->ZOrder(); }
                        pm["z_order"] = picojson::value(sp->ZOrder());
                    }
                    
                    if (sp->Layer() != 0) {
                        if (zoom>=0) { key_strm << "layer=" << sp->Layer(); }
                        pm["layer"] = picojson::value(sp->Layer());
                    }
                    pm["way_area"] = picojson::value(sp->Area());
                }
                    
                std::get<1>(entry)=key_strm.str();
                std::get<3>(entry) = (zoom==-1);
                result.push_back(entry);
            }
                
            return result;
        }
    
        
            
    private:
        //std::shared_ptr<geos_base> gs;
        feature_spec features;
        std::map<std::string,extra_tags_spec> extra_tags;
        int64 max_tile_;
};
        
        
    
    


class MakeTileData {
    public:
        MakeTileData(
            feature_spec features,
            std::map<std::string,extra_tags_spec> extra_tags,
            int64 max_tile, bool clip_geoms_, bool use_obj_tile_,
            bbox filter_box_, bool fix_geoms_, bool simplify_max_)
            
        : featureproperties(features, extra_tags, max_tile), clip_geoms(clip_geoms_), use_obj_tile(use_obj_tile_), filter_box(filter_box_), fix_geoms(fix_geoms_), simplify_max(simplify_max_) {
            
            
            
        }
        
        
        virtual ~MakeTileData() {}
        
        void make_alt_feature_data(xyz basetile, std::shared_ptr<geos_base> gst, std::shared_ptr<oqt::BaseGeometry> geom, const std::function<void(xyz,std::string,alt_feature_data&&)>& cb) const {
        
            oqt::TimeSingle ts;
            int64 minzoom = geom->MinZoom();
            if (minzoom < 0) { throw std::domain_error("not a valid geom"); }
            
            int64 id_ = geom->Id();
            if (geom->Type()==oqt::ElementType::ComplicatedPolygon) {
                id_ *= -1;
            }
            xyz feat_basetile = basetile;
            if (use_obj_tile) {
                auto ob = geom->Bounds();
                int64 qt=oqt::quadtree::calculate(ob.minx,ob.miny,ob.maxx,ob.maxy,0,featureproperties.max_tile());
                
                feat_basetile = oqt::quadtree::tuple(qt);
            }
            
            size_t nt=0; double tss=10;
            int64 npmax = 1ull<< (32 - 14);//featureproperties.max_tile());
            for (int64 geom_ty=0; geom_ty < 3; geom_ty++) {
                
                                    
                            
                
                auto tabs = featureproperties.prep_feature_table(geom, geom_ty);
                
                if (!tabs.empty()) {
                    if ((geom_ty==2) && (geom->Type() == oqt::ElementType::ComplicatedPolygon)) {
                        auto cp = std::dynamic_pointer_cast<oqt::geometry::ComplicatedPolygon>(geom);
                        if (cp->Area() > 200000000000) {
                            oqt::Logger::Message lm;
                            lm << "skip large poly: " << -cp->Id() << " area=" << cp->Area();
                            for (const auto& tg: cp->Tags()) {
                                if ((tg.key == "name") || (tg.key=="place") || (tg.key=="boundary")) {
                                    lm << " " << tg.key << "='" << tg.val << "'";
                                }
                            }
                            continue;
                        }
                    } else if ((geom_ty==2) && (geom->Type() == oqt::ElementType::SimplePolygon)) {
                        auto sp = std::dynamic_pointer_cast<oqt::geometry::SimplePolygon>(geom);
                        if (sp->Area() > 200000000000) {
                            oqt::Logger::Message lm;
                            lm << "skip large poly: " << sp->Id() << " area=" << sp->Area();
                            for (const auto& tg: sp->Tags()) {
                                if ((tg.key == "name") || (tg.key=="place") || (tg.key=="boundary")) {
                                    lm << " " << tg.key << "='" << tg.val << "'";
                                }
                            }
                            continue;
                        }
                    }
                    
                    
                    auto orig = make_geometry(gst, filter_box, geom, geom_ty);
                    if (!orig) { 
                        continue;
                    }
                    auto bnds=orig->bounds();
                    
                    int64 mxt=14;//featureproperties.max_tile();
                    
                    auto makefeat = [&cb, id_, minzoom, &bnds, &tabs, mxt, npmax,&ts,&nt,&geom,&tss](xyz t, std::shared_ptr<geos_geometry> g) {
                        
                        if (ts.since()>tss) {
                            std::cout << "make_alt_feature_data " << geom->Type() << " " << geom->Id() << ": " << ts << " tile " << nt << " " << tup_str(t) << std::endl;
                            tss+=10;
                        }
                                                
                        if (tabs.count(std::get<2>(t))>0) {
                            
                                                        
                            for (const auto& ft: tabs.at(std::get<2>(t))) {
                                
                                auto ss = prep_alt_feature_data(id_, std::get<2>(ft), g, minzoom, std::get<2>(t)==mxt ? npmax : 4096 );
                                cb(t, std::get<0>(ft), std::move(ss));
                                nt++;
                            }
                            
                            //cb(t, alt_feature_data{id_,tabsp.at(std::get<2>(t)),minzoom,g,bnds});
                        }
                    };
                    
                    prep_geometries_clip_cb(orig, feat_basetile, minzoom,featureproperties.max_tile(), fix_geoms,false, makefeat, simplify_max);
                    
                    
                }
            }
            if (ts.since()>10) {
                std::cout << "make_alt_feature_data " << geom->Type() << " " << geom->Id() << ": " << ts << " tile " << nt << " done" << std::endl;
            }
        
        }
        
        std::shared_ptr<alt_tile_data> make_alt_tile_data(oqt::PrimitiveBlockPtr block) const {
            if (!block) { throw std::domain_error("no block"); }
            
            auto gst = std::make_shared<geos_base>();
            
            
            xyz tt = oqt::quadtree::tuple(block->Quadtree());
            
            auto result = std::make_shared<alt_tile_data>();
            result->tile_xyz = tt;
            
            auto cb = [&result](xyz tl, std::string tb, alt_feature_data&& f) {
                result->objs[std::make_pair(tl,tb)].push_back(std::move(f));
            };
            
            for (auto ele: block->Objects()) {
                auto geom = std::dynamic_pointer_cast<oqt::BaseGeometry>(ele);
                if (geom && (geom->MinZoom()>=0)) {
                    
                    try {
                        make_alt_feature_data(tt, gst,geom, cb);
                        
                    } catch (std::exception& ex) {
                        std::cout << "make_alt_feature_data(" << ele->Type() << " " << ele->Id() << " failed" << ex.what() << std::endl;
                    }
                }
            }
            return result;
        }
        
        std::shared_ptr<geos_geometry> simplify_geom(std::shared_ptr<geos_geometry> geom, int64 zoom) const {
            //if (zoom >= featureproperties.max_tile()) { return geom; }
            
            return simplify_geom_check(geom, zoom);
            
        }
        
        std::pair<std::map<xyz,std::shared_ptr<geos_geometry>>,bbox> prep_geometries(
            std::shared_ptr<geos_base> gst, xyz basetile, std::shared_ptr<oqt::BaseGeometry> geom, int64 ty, int64 minzoom, bool simp_max) const {
            
            
            auto orig = make_geometry(gst, filter_box, geom, ty);
            if (!orig) {
                return std::make_pair(std::map<xyz,std::shared_ptr<geos_geometry>>(), bbox());
            }
            if (clip_geoms) {
                return std::make_pair(prep_geometries_clip(orig, basetile, minzoom, featureproperties.max_tile(), fix_geoms,false, simp_max), orig->bounds());
            }
            
            return std::make_pair(prep_geometries_noclip(orig, minzoom, featureproperties.max_tile(),simp_max), orig->bounds());
        }
        
        
        const FeatureProperties& get_feature_properties() const { return featureproperties; }
        
    private:
        FeatureProperties featureproperties;
        bool clip_geoms;
        bool use_obj_tile;
        bbox filter_box;
        bool fix_geoms;
        bool simplify_max;
};
    

/*

class GroupAltTiles {
    public:
        GroupAltTiles(const std::vector<xyz>& qts_in, int64 maxlevel_) : maxlevel(maxlevel_) {
            for (auto qt: qts_in) {
                int64 q=oqt::quadtree::from_tuple(std::get<0>(qt),std::get<1>(qt),std::get<2>(qt));
                qts.insert(oqt::quadtree::round(q, maxlevel));
            }
            
        }
        
        
        virtual ~GroupAltTiles() {}
        
       
        xyz find_key(xyz tl) {
            int64 q=oqt::quadtree::from_tuple(std::get<0>(tl),std::get<1>(tl),std::get<2>(tl));
            auto it = qts.upper_bound(q);
            if (it==qts.begin()) {
                return xyz{0,0,0};
            }
            it--;
            return oqt::quadtree::tuple(*it);
        }
        
        std::shared_ptr<alt_tile_data_vec> call(std::shared_ptr<alt_tile_data> val) {
            if (!val) { return nullptr; }
            std::map<xyz, std::shared_ptr<alt_tile_data>> temp;
            for (const auto& gg: val->objs) {
                auto t = find_key(gg.first.first);
                if (temp.count(t)==0) {
                    temp[t] = std::make_shared<alt_tile_data>();
                    temp[t]->tile_xyz = t;
                }
                temp[t]->objs.insert(gg);
            }
            
            auto res = std::make_shared<alt_tile_data_vec>();
            res->tile_xyz = val->tile_xyz;
            res->tiles.reserve(temp.size());
            for (const auto& tt: temp) {
                res->tiles.push_back(tt.second);
            }
            return res;
        }
    private:
        //std::set<xyz> qts;
        std::set<int64> qts;
        int64 maxlevel;
};
    
*/



std::function<void(oqt::PrimitiveBlockPtr)> make_maketiledata_alt_callback(
feature_spec features, std::map<std::string,extra_tags_spec> extra_tags,
    int64 max_tile, bool use_obj_tile,
    bbox filter_box, bool fix_geoms, bool simplify_max, std::function<void(std::shared_ptr<alt_tile_data>)> cb) {
    
    try {
        auto mtd = std::make_shared<MakeTileData>(features, extra_tags, max_tile, true, use_obj_tile, filter_box, fix_geoms, simplify_max);
        
        return [mtd, cb](oqt::PrimitiveBlockPtr bl) {
            if (!bl) { cb(nullptr); return; }
            
            auto res=mtd->make_alt_tile_data(bl);
            cb(res);
        };
    } catch (std::exception& ex) {
        std::cout << "MakeTileData failed " << ex.what() << std::endl;
    }
    return nullptr;
}

/*
std::function<void(oqt::PrimitiveBlockPtr)> make_maketiledata_groupalt_callback(
feature_spec features, std::map<std::string,extra_tags_spec> extra_tags,
    int64 max_tile, bool use_obj_tile,
    bbox filter_box, bool fix_geoms, bool simplify_max,
    const std::vector<xyz>& qts, std::function<void(std::shared_ptr<alt_tile_data_vec>)> cb) {
    
    try {
        auto mtd = std::make_shared<MakeTileData>(features, extra_tags, max_tile, true, use_obj_tile, filter_box, fix_geoms,simplify_max);
        auto gat = std::make_shared<GroupAltTiles>(qts,max_tile);
        return [mtd, cb, gat](oqt::PrimitiveBlockPtr bl) {
            if (!bl) { cb(nullptr); return; }
            
            auto res=mtd->make_alt_tile_data(bl);
            auto gg = gat->call(res);
            cb(gg);
        };
    } catch (std::exception& ex) {
        std::cout << "MakeTileData failed " << ex.what() << std::endl;
    }
    return nullptr;
}
*/


void export_maketiledata(py::module& m) {
    
    
    
    py::class_<MakeTileData>(m, "MakeTileData")
        .def(py::init<feature_spec,std::map<std::string,extra_tags_spec>,int64,bool,bool,bbox,bool,bool>())
        //.def("make_feature_data", &MakeTileData::make_feature_data)
        //.def("make_tile_data", &MakeTileData::make_tile_data)
        
        .def("prep_geometries", &MakeTileData::prep_geometries)
        
        .def("prep_feature_table", [](const MakeTileData& mtd, std::shared_ptr<oqt::BaseGeometry> geom, int64 geom_ty) {
            return mtd.get_feature_properties().prep_feature_table(geom,geom_ty);
        })
        .def("check_feature", [](const MakeTileData& mtd, std::shared_ptr<oqt::BaseGeometry> geom, int64 zoom) {
            py::gil_scoped_release r;
            
            const auto& fp=mtd.get_feature_properties();
            
            if (!geom) { throw std::domain_error("not a geom??"); }
            
            std::map<std::string,std::vector<oqt::Tag>> tabs;
            int t=geom_type(geom);
            for (const auto& tg: geom->Tags()) {
                fp.check_feature(tabs, t, tg, zoom);
            }
            return tabs;
        })
        
        .def("make_alt_feature_data", &MakeTileData::make_alt_feature_data)
        
        .def("make_alt_tile_data", &MakeTileData::make_alt_tile_data)
    ;
    
    m.def("make_geometry", &make_geometry);
    m.def("prep_geometries_clip", &prep_geometries_clip);
    m.def("prep_geometries_noclip", &prep_geometries_noclip);
    m.def("simplify_geom_check", &simplify_geom_check);
    m.def("clip_geom_tile", &clip_geom_tile);
    
    
    m.def("make_maketiledata_alt_callback", &make_maketiledata_alt_callback);
    
    //m.def("make_maketiledata_groupalt_callback", &make_maketiledata_groupalt_callback);
     
    py::class_<alt_tile_data, std::shared_ptr<alt_tile_data>>(m, "alt_tile_data")
        .def(py::init<>())
        .def_readonly("tile_xyz", &alt_tile_data::tile_xyz)
        .def_readonly("objs", &alt_tile_data::objs)
        .def("add", [](alt_tile_data& atd, xyz tile, std::string tab, alt_feature_data f) {
            atd.objs[std::make_pair(tile,tab)].push_back(f);
        })
    ;
    /*py::class_<GroupAltTiles,std::shared_ptr<GroupAltTiles>>(m,"GroupAltTiles")
        .def(py::init<std::vector<xyz>,int64>())
        .def("__call__", &GroupAltTiles::call)
        .def("find_key", &GroupAltTiles::find_key);
    ;*/
    m.def("is_box", &is_box);
}

