#include "mvt.hpp"
#include "maketiledata.hpp"
#include <oqt/utils/compress.hpp>
#include <oqt/utils/pbf/packedint.hpp>
#include <oqt/utils/pbf/protobuf.hpp>
#include <oqt/utils/pbf/varint.hpp>
#include <oqt/utils/pbf/fixedint.hpp>


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>
namespace py = pybind11;
using namespace std::literals::string_literals;

PYBIND11_DECLARE_HOLDER_TYPE(XX, std::shared_ptr<XX>);


oqt::int64 asint(double v) {
    if (v>=0) {
        return v+0.5;
    }
    return v-0.5;
}    


struct deltazz {
    deltazz(int64 np) : x(0), y(0) {
        sf = np / 256.0;
    }
    int64 x; int64 y;
    double sf;
    size_t operator()(std::vector<oqt::uint64>& result, size_t i, std::pair<double,double> n) {
        int64 xi = asint(sf*n.first);
        int64 yi = asint(sf*(256-n.second));        
        result[i] = oqt::zig_zag(xi-x);
        result[i+1] = oqt::zig_zag(yi-y);
        x=xi;
        y=yi;
        return i+2;
        
    }
};
class dedeltazz {
    public:
        dedeltazz(int64 np) : x(0), y(0) {
            sf = 256.0 / np;
        }
        
        std::pair<double,double> operator()(const std::vector<oqt::uint64>& cmds, size_t pos) {
            x+=oqt::un_zig_zag(cmds[pos]);
            y+=oqt::un_zig_zag(cmds[pos+1]);
            return std::make_pair(x*sf, 256-y*sf);
        }

    private:
        int64 x,y;
        double sf;
};

std::string write_points(const std::vector<std::shared_ptr<geos_geometry>>& pts, int64 np) {
    std::vector<oqt::uint64> result(1+pts.size()*2);
    result[0] = 1 | (pts.size()<<3);
    deltazz zz(np);
    size_t i=1;
    
    for (const auto& g: pts) {
        auto cc = g->get_coords();
        if (cc.size()!=1) { throw std::domain_error("??"); }
        i = zz(result, i, cc[0]);
    }
    return oqt::write_packed_int(result);
}

size_t write_coords(std::vector<oqt::uint64>& result, size_t pos, deltazz& zz, const coords& cc, bool closed) {
    if (cc.size() < 2) { throw std::domain_error("??"); }
    if (closed && (cc.size() < 3) && (cc.front()!=cc.back())) { throw std::domain_error("??"); }
        
    result[pos] = 1 | (1<<3);
    pos++;
    
    pos = zz(result, pos, cc[0]);
    
    size_t npts = cc.size()-1 - (closed? 1: 0);
    result[pos] = 2 | (npts<<3);
    pos++;
    
    for (size_t i=0; i < npts; i++) {
        pos = zz(result, pos, cc[1+i]);
    }
    
    if (closed) {
        result[pos] = 7 | (1<<3);
        pos++;
    }
    return pos;
}


size_t coords_len(const coords& cc, bool closed) {
    if (closed) {
        return 1 + 2 + 1 + 2*(cc.size()-2) + 1;
    }
    
    return 1 + 2 + 1 + 2*(cc.size()-1);
}

std::string write_linestrings(const std::vector<std::shared_ptr<geos_geometry>> lines, int np) {
    deltazz zz(np);
    std::vector<oqt::uint64> result;
    size_t pos=0;
    for (const auto& l: lines) {
        auto cc = l->get_coords();
        result.resize(pos + coords_len(cc,false));
        pos=write_coords(result, pos, zz, cc, false);
    }
    
    return oqt::write_packed_int(result);
}
        
std::string write_polygons(const std::vector<std::shared_ptr<geos_geometry>> polygons, int np) {
    deltazz zz(np);
    std::vector<oqt::uint64> result;
    size_t pos=0;
    for (const auto& p: polygons) {
        auto rings = p->get_rings();
        size_t ll = 0;
        for (const auto& cc : rings) { ll += coords_len(cc, true); }
        result.resize(pos + ll);
        for (const auto& cc: rings) {
            pos=write_coords(result, pos, zz, cc, true);
        }
    }
    
    return oqt::write_packed_int(result);
}

std::pair<size_t,std::string> pack_geometry(std::shared_ptr<geos_geometry> geom, int64 np) {
    if ((geom->type_id()== 0) || (geom->type_id()==4) ) {
        return std::make_pair(1, write_points(geom->get_geometries(), np));
    }
    
    if ((geom->type_id()== 1) || (geom->type_id()==5) ) {
        return std::make_pair(2, write_linestrings(geom->get_geometries(), np));
    }
    
    if ((geom->type_id()== 3) || (geom->type_id()==6) ) {
        return std::make_pair(3, write_polygons(geom->get_geometries(), np));
    }
    
    return std::make_pair(4, geom->wkb());
    
}

bool write_geometry(std::list<oqt::PbfTag>& msgs, std::shared_ptr<geos_geometry> geom, int64 np) {
    
    size_t gt; std::string gd;
    std::tie(gt,gd) = pack_geometry(geom, np);
    
    msgs.push_back(oqt::PbfTag{3, gt, ""});
    msgs.push_back(oqt::PbfTag{4, 0, gd});
    return (gt==1) || (gt==2) || (gt==3);
    
}


coords read_coords(const std::vector<oqt::uint64>& src, size_t& pos, dedeltazz& zz, bool close) {
    if (pos >= src.size()) { throw std::domain_error("out of range"); }
    if (src[pos] != (1 | (1<<3))) { throw std::domain_error("ring doesn't start with a single moveto"); }
    
    if ((src[pos+3] & 7) != 2) { throw std::domain_error("ring doesn't have a lineto"); }
    size_t nlt = src[pos+3]>>3;
    
    if (close) {
        if ((src.size() < (pos+4+2*nlt+1)) || (src[pos+4+2*nlt] != 15)) {
            throw std::domain_error("close cmd not present");
        }
    }
    
    coords result(1+nlt+(close?1:0));
    
    pos++;
    result[0] = zz(src, pos);
    pos+=3;
    for (size_t i=0; i < nlt; i++) {
        result[i+1] = zz(src, pos);
        pos+=2;
    }
    if (close) {
        result[1+nlt] = result[0];
        pos++;
    }
    
    return result;
}

coords read_coords_points(const std::vector<oqt::uint64>& cmds, size_t& pos, dedeltazz& zz) {
    if (cmds.size()<(3+pos)) { throw std::domain_error("no points??"); }
    if ((cmds[pos] & 7)!=1) { throw std::domain_error("points don't start with a move to??"); }
    
    coords geoms;
    
    size_t num_points = cmds[pos]>>3;
    geoms.reserve(num_points);
    
    pos++;
    
    
    for (size_t i=0; i < num_points; i++) {
        geoms.push_back(zz(cmds, pos));
        pos+=2;
        
    }
    return geoms;
}

std::shared_ptr<geos_geometry> read_points(std::shared_ptr<geos_base> gs, const std::string& data, int64 np) {
    
    std::vector<oqt::uint64> cmds = oqt::read_packed_int(data);
    std::vector<std::shared_ptr<geos_geometry>> geoms;
    dedeltazz zz(np);
    size_t pos=0;
    while (pos < cmds.size()) {
        auto cc = read_coords_points(cmds, pos, zz);
        geoms.reserve(geoms.size()+cc.size());
        for (const auto& p: cc) {
            geoms.push_back(geos_geometry_create_point(gs, p.first, p.second));
        }
    }
    
    return geos_geometry_create_collection(gs, geoms);
}
    
std::shared_ptr<geos_geometry> read_linestrings(std::shared_ptr<geos_base> gs, const std::string& data, int64 np) {
    
    std::vector<oqt::uint64> cmds = oqt::read_packed_int(data);
    std::vector<std::shared_ptr<geos_geometry>> geoms;
    dedeltazz zz(np);
    
    size_t pos=0;
    while (pos < cmds.size()) {
        coords ln = read_coords(cmds, pos, zz,false);
        geoms.push_back(geos_geometry_create_linestring(gs, ln));
    }
    
    return geos_geometry_create_collection(gs, geoms);
}


double ring_area(const coords& r) {
    if (r.size()<3) { throw std::domain_error("not a ring"); }
    
    
    double a=0;
    for (size_t i=0; i < (r.size()-1); i++) {
        a+=(r[i+1].first-r[i].first)*(r[i+1].second+r[i].second);
    }
    return a/2;
}

std::shared_ptr<geos_geometry> read_polygons(std::shared_ptr<geos_base> gs, const std::string& data, int64 np) {
    
    std::vector<oqt::uint64> cmds = oqt::read_packed_int(data);
    std::vector<std::shared_ptr<geos_geometry>> geoms;
    dedeltazz zz(np);
    
    std::vector<coords> rings;
    
    std::vector<std::shared_ptr<geos_geometry>> polys;
    
    size_t pos=0;
    while (pos < cmds.size()) {
        coords r = read_coords(cmds, pos, zz,true);
        
        if (ring_area(r)>0) {
            if (!rings.empty()) {
                polys.push_back(geos_geometry_create_polygon(gs, rings));
            }
            rings.clear();
        }
        rings.push_back(r);
        //geoms.push_back(geos_geometry_create_linestring(gs, ln));
    }
    polys.push_back(geos_geometry_create_polygon(gs, rings));
    
    if (polys.size()==1) {
        return polys[0];
    }
    auto p= geos_geometry_create_collection(gs, polys);
    if (!p->is_valid()) {
        return p->unary_union();
    }
    return p;
    
}

std::shared_ptr<geos_geometry> read_geometry(std::shared_ptr<geos_base> gs, oqt::uint64 gt, const std::string& data, int64 np) {
    
    if (gt==1) { return read_points(gs, data, np); }
    if (gt==2) { return read_linestrings(gs, data, np); }
    if (gt==3) { return read_polygons(gs, data, np); }
    if (gt==4) { return geos_geometry_create_wkb(gs, data); }
    
    throw std::domain_error("unexpected type "+std::to_string(gt));
    return nullptr;
}
    
    
typedef std::map<std::string,oqt::uint64> keylookup;


oqt::uint64 find_key(std::list<oqt::PbfTag>& layermsgs, size_t tg, keylookup& keys, const std::string& s) {
    
    auto it = keys.find(s);
    if (it!=keys.end()) { return it->second; }
    
    layermsgs.push_back(oqt::PbfTag{tg, 0, s});
    
    oqt::uint64 l = keys.size();
    keys.insert(std::make_pair(s, l));
    
    return l;
}

std::string pack_value(const picojson::value& v) {
    
    if (v.is<std::string>()) {
        return oqt::pack_pbf_tags({oqt::PbfTag{1,0,v.get<std::string>()}});
    }
    
    if (v.is<int64>()) {
        return oqt::pack_pbf_tags({oqt::PbfTag{6,oqt::zig_zag(v.get<int64>()),""}});
    }
    if (v.is<double>()) {
        
        std::string r(9,'\0');
        r[0] = (3<<3) | 1;
        oqt::write_double_le(r, 1, v.get<double>());
        return r;
    }
    if (v.is<bool>() && v.get<bool>()) {
        return oqt::pack_pbf_tags({oqt::PbfTag{7,1,""}});
    }
    return oqt::pack_pbf_tags({oqt::PbfTag{7,0,""}});
}



double read_double(const std::string& data, size_t pos) {
    union {
        double d;
        char sc[8];
    };
    
    sc[0] = data[pos];
    sc[1] = data[pos+1];
    sc[2] = data[pos+2];
    sc[3] = data[pos+3];
    sc[4] = data[pos+4];
    sc[5] = data[pos+5];
    sc[6] = data[pos+6];
    sc[7] = data[pos+7];
    return d;
}

double read_float(const std::string& data, size_t pos) {
    union {
        float f;
        char sc[4];
    };
    
    sc[0]=data[pos];
    sc[1]=data[pos+1];
    sc[2]=data[pos+2];
    sc[3]=data[pos+3];
    
    return f;
}

picojson::value read_value(const std::string& data) {
    
    //tg.tag==2: float 32bit
    if ( (data.size()==5) && (data[0] == ((2<<3)|5))) { 
        double r=read_float(data,1);
        return picojson::value(r);
        throw std::domain_error("can't read value "+data);
    } 
    
    //tg.tag==3: double 64bit
    if ( (data.size()==9) && (data[0] == ((3<<3)|1))) {
        double r=read_double(data,1);
        return picojson::value(r);
    } 
    size_t pos=0;
    oqt::PbfTag tg = oqt::read_pbf_tag(data,pos);
    if (tg.tag==1) { return picojson::value(tg.data); }
    if (tg.tag==4) { return picojson::value(int64(tg.value)); }
    if (tg.tag==5) { return picojson::value(int64(tg.value)); }
    if (tg.tag==6) { return picojson::value(oqt::un_zig_zag(tg.value)); }
    if (tg.tag==7) { return picojson::value(tg.value==1); }
    
    throw std::domain_error("can't read value "+data);
}

py::object read_value_python(const std::string& data) {
    //tg.tag==2: float 32bit
    if ( (data.size()==5) && (data[0] == ((2<<3)|5))) { 
        double r=read_float(data,1);
        return py::cast(r);
    } 
    
    //tg.tag==3: double 64bit
    if ( (data.size()==9) && (data[0] == ((3<<3)|1))) {
        double r=read_double(data,1);
        return py::cast(r);
    } 
    size_t pos=0;
    oqt::PbfTag tg = oqt::read_pbf_tag(data,pos);
    if (tg.tag==1) { return py::cast(tg.data); }
    if (tg.tag==4) { return py::cast(tg.value); }
    if (tg.tag==5) { return py::cast(int64(tg.value)); }
    if (tg.tag==6) { return py::cast(oqt::un_zig_zag(tg.value)); }
    if (tg.tag==7) { return py::cast(tg.value==1); }
    
    throw std::domain_error("can't read value "+data);
}

double read_value_double(const std::string& data) {
    //tg.tag==3: double 64bit
    if ( (data.size()==9) && (data[0] == ((3<<3)|1))) {
        return read_double(data,1);
        
    } 
    size_t pos=0;
    oqt::PbfTag tg = oqt::read_pbf_tag(data,pos);
    
    if (tg.tag==4) { return (double) tg.value; }
    if (tg.tag==5) { return (double) tg.value; }
    if (tg.tag==6) { return (double) oqt::un_zig_zag(tg.value); }
    
    throw std::domain_error("can't read value "+data+" as double");
}
int64 read_value_integer(const std::string& data) {
    size_t pos=0;
    oqt::PbfTag tg = oqt::read_pbf_tag(data,pos);

    if (tg.tag==4) { return int64(tg.value); }
    if (tg.tag==5) { return int64(tg.value); }
    if (tg.tag==6) { return oqt::un_zig_zag(tg.value); }
    throw std::domain_error("can't read value "+data+" as integer");
}


oqt::uint64 find_val(std::list<oqt::PbfTag>& layermsgs, oqt::uint64 tg, keylookup& vals, const picojson::value& v) {
    
    std::string s = pack_value(v);
    return find_key(layermsgs, tg, vals, s);
}

std::string tile_fn(int64 x, int64 y, int64 z) {
    std::stringstream ss;
    ss << x << "_" << y << "_" << z << ".mvt";
    return ss.str();
}

std::string pack_property_interim(const std::string& k, const picojson::value& v) {
    std::list<oqt::PbfTag> msgs;
    msgs.push_back(oqt::PbfTag{1, 0, k});
    msgs.push_back(oqt::PbfTag{2, 0, pack_value(v)});
    return oqt::pack_pbf_tags(msgs);
}

    
std::pair<std::string,double> pack_propertymap_interim(const property_map& pm, oqt::uint64 gt) {
    std::list<oqt::PbfTag> msgs;
    double areaorlen=0;
    for (const auto& kv: pm) {
        if ((gt==2) && (kv.first=="way_length"s) && (kv.second.is<double>())) {
            areaorlen = kv.second.get<double>();
        } else if ((gt==3) && (kv.first=="way_area"s) && (kv.second.is<double>())) {
            areaorlen = kv.second.get<double>();
        } else {
            msgs.push_back(oqt::PbfTag{1, 0, pack_property_interim(kv.first, kv.second)});
        }
    }
    return std::make_pair(oqt::pack_pbf_tags(msgs),areaorlen);
}


alt_feature_data prep_alt_feature_data(int64 id, const property_map& pm, std::shared_ptr<geos_geometry> geom, int64 minzoom, int64 np) {
    
    alt_feature_data result;
    result.id=id;
    result.minzoom=minzoom;
    result.np = np;
    
    
    std::tie(result.geom_type, result.geom_data) = pack_geometry(geom, np);
    std::tie(result.properties,result.areaorlen) = pack_propertymap_interim(pm, result.geom_type);
    
    return result;
}


void pack_feature_from_interim(std::list<oqt::PbfTag>& layermsgs, keylookup& keys, keylookup& vals, const alt_feature_data& f, int64 np) {

    std::list<oqt::PbfTag> featmsgs;
    featmsgs.push_back(oqt::PbfTag{1, (oqt::uint64) f.id, ""});
    
    
    
    std::vector<oqt::uint64> kvs;
    
    size_t pos=0;
    oqt::PbfTag tg = oqt::read_pbf_tag(f.properties, pos);
    for ( ; tg.tag>0; tg = oqt::read_pbf_tag(f.properties, pos)) {
        if (tg.tag!=1) { throw std::domain_error("???"); }
        size_t pos2=0;
        auto key = oqt::read_pbf_tag(tg.data, pos2).data;
        kvs.push_back(find_key(layermsgs, 3, keys, key));
        auto val = oqt::read_pbf_tag(tg.data, pos2).data;
        kvs.push_back(find_key(layermsgs, 4, vals, val));
        
    }   
    
    
    if ((f.geom_type==2) && (f.areaorlen!=0.0)) {
        kvs.push_back(find_key(layermsgs, 3, keys, "way_length"s));
        kvs.push_back(find_key(layermsgs, 4, vals, pack_value(picojson::value(f.areaorlen))));
    }
    if ((f.geom_type==3) && (f.areaorlen!=0.0)) {
        kvs.push_back(find_key(layermsgs, 3, keys, "way_area"s));
        kvs.push_back(find_key(layermsgs, 4, vals, pack_value(picojson::value(f.areaorlen))));
    }
    if (f.minzoom>=0) {
        kvs.push_back(find_key(layermsgs, 3, keys, "min_zoom"s));
        kvs.push_back(find_key(layermsgs, 4, vals, pack_value(picojson::value(f.minzoom))));
    }
    featmsgs.push_back(oqt::PbfTag{2, 0, oqt::write_packed_int(kvs)});
    featmsgs.push_back(oqt::PbfTag{3, f.geom_type, ""});
    if (f.np!=np) {
        //std::string gd2 = fix_geometry_np(f.geom_type, f.geom_data, f.np, np);
        //featmsgs.push_back(oqt::PbfTag{4, 0, gd2});
        throw std::domain_error("not implemented");
    } else {
        featmsgs.push_back(oqt::PbfTag{4, 0, f.geom_data});
    }
    
    layermsgs.push_back(oqt::PbfTag{2, 0, oqt::pack_pbf_tags(featmsgs)});
    
}

std::string pack_layer_from_interim(int64 tx, int64 ty, int64 tz, const std::string& table, const std::vector<alt_feature_data>& feats, int64 np) {
    
    keylookup keys;
    keylookup vals;

    std::list<oqt::PbfTag> layermsgs;
    layermsgs.push_back(oqt::PbfTag{15,2,""});
    layermsgs.push_back(oqt::PbfTag{1,0,table});
    layermsgs.push_back(oqt::PbfTag{5,(oqt::uint64) np,""});
    
    
    
    for (const auto& f: feats) {
        try {
            pack_feature_from_interim(layermsgs, keys, vals, f, np);
        } catch (std::exception& ex) {
            std::cout << "pack_feature_from_interim " << tup_str(xyz{tx,ty,tz}) << " " << table << " " << " " << f.id << " " << f.properties << " failed: " << ex.what() << std::endl;
        }
    }
    
    return oqt::pack_pbf_tags(layermsgs);
}
std::string pack_tile_mvt_from_interim(int64 x, int64 y, int64 z, const std::map<std::string,std::vector<alt_feature_data>>& tl, int64 np) {
    
    //PackGeom pg(x,y,z,np);
    
    std::list<oqt::PbfTag> msgs;
    for (const auto& st: tl) {
        msgs.push_back(oqt::PbfTag{3, 0, pack_layer_from_interim(x,y,z,st.first, st.second, np)});
    }
    
    return oqt::pack_pbf_tags(msgs);
}

std::shared_ptr<packed_tiles> pack_tiles_mvt_from_interim(const std::map<xyz,std::map<std::string,std::vector<alt_feature_data>>>& tabs, bool compress, int64 maxlevel) {

    auto rr=std::make_shared<packed_tiles>();
    for (const auto& pp: tabs) {
        if (!pp.second.empty()) {
            int64 x=std::get<0>(pp.first); int64 y=std::get<1>(pp.first); int64 z=std::get<2>(pp.first);
            
            std::string packed = pack_tile_mvt_from_interim(x, y, z, pp.second, (z==maxlevel) ? (1ll<<(32-z)) : 4096);
            if (compress) {
                packed = oqt::compress_gzip(tile_fn(x,y,z), packed);
            }
            rr->insert(std::make_pair(std::make_tuple(x,y,z), packed));
        }
    }
    
    return rr;
}

std::string merge_points(const std::vector<std::string>& geoms, int64 np) {
    
    
    
    coords points;
    points.reserve(geoms.size());   
    
    for (const auto& g: geoms) {
        size_t pos=0;
        auto src = oqt::read_packed_int(g);
        
        dedeltazz dzz(np);
        while (pos < src.size()) {
            auto cc = read_coords_points(src, pos, dzz);
            std::copy(cc.begin(), cc.end(), std::back_inserter(points));
        }
        
    }
    
    deltazz zz(np);
    std::vector<oqt::uint64> cmds;
    cmds.resize(1+2*points.size());
    cmds[0] = 1 | (points.size() << 3);
    size_t pos=1;
    for (auto p: points) {
        pos=zz(cmds, pos, p);
    }
    
    return oqt::write_packed_int(cmds);
}

std::string merge_linestrings(const std::vector<std::string>& geoms, int64 np) {
    
    std::vector<coords> ll;
    
    for (const auto& g: geoms) {
        size_t pos=0;
        auto src = oqt::read_packed_int(g);
        
        dedeltazz dzz(np);
        while(pos < src.size()) {
            auto cc = read_coords(src, pos, dzz,false);
            ll.push_back(cc);
        }
    }
    
    deltazz zz(np);
    std::vector<oqt::uint64> result;
    size_t pos=0;
    for (const auto& cc: ll) {
        result.resize(pos + coords_len(cc,false));
        pos=write_coords(result, pos, zz, cc, false);
    }
    
    return oqt::write_packed_int(result);
}

std::string merge_polygons(const std::vector<std::string>& geoms, int64 np) {
    
    std::vector<coords> ll;
    
    for (const auto& g: geoms) {
        size_t pos=0;
        auto src = oqt::read_packed_int(g);
        
        dedeltazz dzz(np);
        while(pos < src.size()) {
            auto cc = read_coords(src, pos, dzz,true);
            ll.push_back(cc);
        }
    }
    
    deltazz zz(np);
    std::vector<oqt::uint64> result;
    size_t pos=0;
    for (const auto& cc: ll) {
        result.resize(pos + coords_len(cc,true));
        pos=write_coords(result, pos, zz, cc, true);
    }
    
    return oqt::write_packed_int(result);
}
    
        
std::string merge_geometries(size_t gt, const std::vector<std::string>& geoms, int64 np) {
    
    if (gt==1) { return merge_points(geoms, np); }
    if (gt==2) { return merge_linestrings(geoms, np); }
    if (gt==3) { return merge_polygons(geoms, np); }
    throw std::domain_error("unexpceted geom type");
    return "";
}

std::vector<alt_feature_data> merge_interim_features(const std::vector<alt_feature_data>& feats, int64 np) {
    
    std::map<std::string,std::pair<alt_feature_data,std::vector<std::string>>> collected;
    
    int64 ii=1;
    
    for (const auto& f: feats) {
        if (f.np!=np) { throw std::domain_error("inconsistent np"); }
        
        std::string pp =f.properties;
        if (f.id<0) {
            pp += ";minus";
        }
        auto it = collected.find(pp);
        
        
        if (it==collected.end()) {
            alt_feature_data nf;
            nf.id = ii++;
            if (f.id<0) { nf.id *= -1; }
            nf.properties=f.properties;
            nf.np = np;
            nf.geom_type= f.geom_type;
            nf.minzoom = f.minzoom;
            nf.areaorlen = f.areaorlen;
            it = collected.insert(std::make_pair(pp, std::make_pair(nf, std::vector<std::string>{f.geom_data}))).first;
        } else {
            
            if (it->second.first.geom_type != f.geom_type) {
                throw std::domain_error("inconsitent geom_type");
            }
            if (f.minzoom < it->second.first.minzoom) {
                it->second.first.minzoom = f.minzoom;
            }
            if (f.areaorlen != 0.0) {
                it->second.first.areaorlen += f.areaorlen;
            }
            
            it->second.second.push_back(f.geom_data);
        }
        
    }
    
    std::vector<alt_feature_data> result;
    result.reserve(collected.size());
    for (auto& fp: collected) {
        if (fp.second.second.size()==0) {
            //???
        } else if (fp.second.second.size()==1) {
            fp.second.first.geom_data=fp.second.second.front();
        } else {
            fp.second.first.geom_data = merge_geometries(fp.second.first.geom_type, fp.second.second, np);
        }
        result.push_back(fp.second.first);
    }
    //std::cout << "merge_interim_features: " << feats.size() << " => " << result.size() << std::endl;
    return result;
}
            
            
inline bool afd_cmp(const alt_feature_data& l, const alt_feature_data& r) {
    if (l.id==r.id) {
        if (l.minzoom==r.minzoom) {
            if (l.properties==r.properties) {
                if (l.areaorlen==r.areaorlen) {
                    return l.geom_data < r.geom_data;
                }
                return l.areaorlen < r.areaorlen;
            }
            return l.properties < r.properties;
        }
        return l.minzoom < r.minzoom;
    }
    return l.id < r.id;
}



std::map<xyz,std::map<std::string,std::vector<alt_feature_data>>> merge_interim_tiles(std::shared_ptr<alt_tile_data_vec> tv, int64 maxlevel, bool sortobjs, bool mergefeats) {
    
    std::map<xyz,std::map<std::string,std::vector<alt_feature_data>>> tabs;
    
    if (!tv) { return tabs; }
    if (tv->tiles.empty()) { return tabs; }
    if (tv->to_finish.empty()) { return tabs; }
    
    for (auto f: tv->to_finish) {
        tabs[f]=std::map<std::string,std::vector<alt_feature_data>>();
    }
    
    for (auto t: tv->tiles) {
        if (!t) { continue; }
        for (const auto& ff: t->objs) {
            if (tabs.count(ff.first.first)==0) {
                //pass
            } else {
                auto& vv = tabs.at(ff.first.first)[ff.first.second];
                vv.reserve(vv.size()+ff.second.size());
                std::copy(ff.second.begin(),ff.second.end(), std::back_inserter(vv));
            }
        }
    }
    
    for (auto& ff: tabs) {
        
        for (auto& st: ff.second) {
            if (sortobjs) {
                std::sort(st.second.begin(), st.second.end(), afd_cmp);
            }
            
        
            if ((std::get<2>(ff.first)<14) && mergefeats) {// && (st.first=="highway"s)) {
                try {
                    auto pp = merge_interim_features(st.second, 4096);
                    st.second.swap(pp);
                } catch (std::exception& ex) {
                    std::cout << "merge_interim_features " << tup_str(ff.first) << " [highway] failed " << ex.what() << std::endl;
                }
                    
                    
            }
        }
    }
    return tabs;
    
}

alt_tile_data_vec_cb make_merge_interim_tiles_callback(bool compress, int64 maxlevel, bool sortobjs, bool mergefeats, std::function<void(std::shared_ptr<packed_tiles>)> cb) {
    
    return [compress, maxlevel, sortobjs, mergefeats, cb](std::shared_ptr<alt_tile_data_vec> tv) {
        if (!tv) {
            cb(nullptr);
            return;
        }
        auto tabs = merge_interim_tiles(tv, maxlevel, sortobjs, mergefeats);
        auto pp = pack_tiles_mvt_from_interim(tabs, compress, 14);
        cb(pp);
    };
}
    
    


py::object cast_packed_tiles(std::shared_ptr<packed_tiles> pt) {
    if (!pt) { return py::none(); }
    py::dict result;
    for (const auto& kv: *pt) {
        result[py::cast(kv.first)]=py::bytes(kv.second);
    }
    return result;
}



void read_mvt_feature(std::vector<alt_feature_data>& features,
    const std::vector<std::string>& keys,
    const std::vector<std::string>& vals,
    size_t np,
    const std::string& data) {
    
    std::string props;
    
    alt_feature_data feat;
    feat.areaorlen=0.;
    feat.np = np;
    
    size_t pos=0;
    auto tg = oqt::read_pbf_tag(data,pos);
    for ( ; tg.tag>0; tg=oqt::read_pbf_tag(data,pos)) {
        if (tg.tag==1) {
            feat.id = (int64) tg.value;
        } else if (tg.tag==2) {
            props = tg.data;
            
            
        } else if (tg.tag==3) {
            feat.geom_type = tg.value;
        } else if (tg.tag==4) {
            feat.geom_data = tg.data;
        }
    }
    
    
    std::list<oqt::PbfTag> pm;
            
    auto kvs = oqt::read_packed_int(props);
    if ((kvs.size()%2) != 0) { throw std::domain_error("unexpected keyvals length"); }
    for (size_t i=0; i < kvs.size(); i+=2) {
        const auto& key = keys.at(kvs[i]);
        const auto& val = vals.at(kvs[i+1]);
        if (key=="min_zoom"s) {
            feat.minzoom = read_value_integer(val);
        } else if ((key=="way_length") && (feat.geom_type==2)) {
            feat.areaorlen=read_value_double(val);
        } else if ((key=="way_area") && (feat.geom_type==3)) {
            feat.areaorlen=read_value_double(val);
        } else {
            pm.push_back(oqt::PbfTag{1,0,oqt::pack_pbf_tags({oqt::PbfTag{1,0,key},oqt::PbfTag{2,0,val}})});
            
        }
    }
    feat.properties = oqt::pack_pbf_tags(pm);            
    
    features.push_back(feat);
    
    
}

void read_mvt_layer(std::map<std::string,std::vector<alt_feature_data>>& layers, const std::string& data) {
    
    std::vector<std::string> keys;
    std::vector<std::string> vals;
    
    
    std::vector<alt_feature_data> features;
    
    size_t pos = 0;
    auto tg = oqt::read_pbf_tag(data,pos);
    
    if (tg.tag!=15) { throw std::domain_error("expected tag 15 first..."); }
    if (tg.value!=2) { throw std::domain_error("incorrect version"); }
    
    tg = oqt::read_pbf_tag(data,pos);
    if (tg.tag!=1) { throw std::domain_error("expected tag 1 second..."); }
    std::string layer_name = tg.data;
    
    tg = oqt::read_pbf_tag(data,pos);
    if (tg.tag!=5) { throw std::domain_error("expected tag 5 third..."); }
    size_t np = tg.value;
    
    tg = oqt::read_pbf_tag(data,pos);
    for ( ; tg.tag > 0; tg=oqt::read_pbf_tag(data,pos)) {
        if (tg.tag == 3) {
            keys.push_back(tg.data);
        } else if (tg.tag==4) {
            vals.push_back(tg.data);
        } else if (tg.tag==2) {
            read_mvt_feature(features, keys, vals, np, tg.data);
        }
    }
    layers.emplace(layer_name, features);
}
            
    
    


std::map<std::string,std::vector<alt_feature_data>> read_mvt_tile(const std::string& data, bool compressed) {
    if (compressed) {
        return read_mvt_tile(oqt::decompress_gzip(data), false);
    }
    
    std::map<std::string,std::vector<alt_feature_data>> layers;
    
    size_t pos = 0;
    auto tg = oqt::read_pbf_tag(data,pos);
    for ( ; tg.tag>0; tg = oqt::read_pbf_tag(data,pos)) {
        if (tg.tag == 3) {
            read_mvt_layer(layers, tg.data);
        }
    }
    return layers;
}


void read_property(py::dict& result, const std::string& data) {
    
    size_t pos = 0;
    auto key = oqt::read_pbf_tag(data,pos);
    auto val = oqt::read_pbf_tag(data,pos);
    
    if (key.tag!=1) { throw std::domain_error("??"); }
    if (val.tag!=2) { throw std::domain_error("??"); }
    
    result[py::cast(key.data)] = read_value_python(val.data);
}
    

py::dict read_properties_packed(const std::string& props) {
    py::dict result;
    
    size_t pos=0;
    auto tg = oqt::read_pbf_tag(props, pos);
    for ( ; tg.tag>0; tg = oqt::read_pbf_tag(props, pos)) {
        if (tg.tag==1) {
            read_property(result, tg.data);
        }
    }
    return result;
}



void export_mvt(py::module& m) {
   
    m.def("pack_geometry", [](std::shared_ptr<geos_geometry> g, int64 np) {
        auto gg = pack_geometry(g,np);
        return py::make_tuple(gg.first, py::bytes(gg.second));
    });
           
    m.def("read_geometry", &read_geometry);
    
    m.def("pack_propertymap_interim", [](const property_map& pm, size_t gt) { auto r=pack_propertymap_interim(pm,gt); return py::make_tuple(py::bytes(r.first),py::cast(r.second));});
    
    m.def("prep_alt_feature_data", &prep_alt_feature_data);
    
    py::class_<alt_feature_data>(m, "alt_feature_data")
        .def_readonly("id", &alt_feature_data::id)
        .def_readonly("geom_type", &alt_feature_data::geom_type)
        .def_readonly("minzoom", &alt_feature_data::minzoom)
        .def_readonly("np", &alt_feature_data::np)
        .def_property_readonly("geom_data", [](const alt_feature_data& d) { return py::bytes(d.geom_data); })
        .def_property_readonly("properties", [](const alt_feature_data& d) { return py::bytes(d.properties); })
        .def_readonly("areaorlen", &alt_feature_data::areaorlen)
        .def_property_readonly("tuple", [](const alt_feature_data& d) {
            return py::make_tuple(d.id,py::bytes(d.properties), d.areaorlen, d.geom_type,d.np, py::bytes(d.geom_data), d.minzoom);
        })
    ;
    m.def("merge_interim_features", &merge_interim_features);
    m.def("pack_layer_from_interim", [](int64 x, int64 y, int64 z, std::string t, const std::vector<alt_feature_data>& ff, int64 np) { return py::bytes(pack_layer_from_interim(x,y,z,t, ff,np)); });
    m.def("pack_tile_mvt_from_interim", [](int64 x, int64 y, int64 z, const std::map<std::string,std::vector<alt_feature_data>>& ff, int64 np) { return py::bytes(pack_tile_mvt_from_interim(x,y,z,ff,np)); });
    m.def("pack_tiles_mvt_from_interim", [](const std::map<xyz,std::map<std::string,std::vector<alt_feature_data>>>& tabs, bool compress, int64 maxlevel) { return cast_packed_tiles(pack_tiles_mvt_from_interim(tabs,compress,maxlevel)); });
    
    m.def("merge_interim_tiles", &merge_interim_tiles);
    
    
    m.def("read_mvt_tile", &read_mvt_tile);
    
    m.def("read_properties", &read_properties_packed);
    
};

