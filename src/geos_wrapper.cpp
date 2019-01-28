#include "oqttiles.hpp"
#include "geos_wrapper.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

namespace py = pybind11;
using namespace std::literals::string_literals;

PYBIND11_DECLARE_HOLDER_TYPE(XX, std::shared_ptr<XX>);



bool bbox_overlaps(const bbox& l, const bbox& r) {
    if (std::get<0>(l) > std::get<2>(r)) { return false; }
    if (std::get<1>(l) > std::get<3>(r)) { return false; }
    if (std::get<0>(r) > std::get<2>(l)) { return false; }
    if (std::get<1>(r) > std::get<3>(l)) { return false; }
    return true;
}

const double earth_width = 40075016.68557849;
bbox tile_bound(int64 x, int64 y, int64 z, double buf) {
    if ((z<0) || (z>20)) { throw std::domain_error("out of range"); }
    
    int64 tz = 1ll << z;
    
    if ((x<0) || (x>=tz)) { throw std::domain_error("out of range"); }
    if ((y<0) || (y>=tz)) { throw std::domain_error("out of range"); }
    
    double tw = earth_width / (1ll<<z);
    
    double x0 = x*tw - earth_width/2;
    double y0 = earth_width/2 - y*tw;
    
    double bv = tw*buf/2;
    
    return std::make_tuple(x0-bv, y0 - tw -bv, x0+tw+bv, y0+bv);
}
double zoom_scale(int64 z) {
    if ((z<0) || (z>20)) { throw std::domain_error("out of range"); }
    return earth_width / (1ll<<z) / 256.0;
}
    


geos_base::geos_base() {
    handle_ = GEOS_init_r();
    GEOSContext_setNoticeMessageHandler_r(handle_, this->notice, data);
    GEOSContext_setErrorMessageHandler_r(handle_, this->error, data);
    lasterror="";
}

geos_base::~geos_base() {
    //std::cout << "~geos_base()" << std::endl;
    GEOS_finish_r(handle_);
}

void geos_base::error(const char* msg, void* data) {
    std::cout << "geos error " << msg << std::endl;
    //lasterror=std::string(msg);
    
}
void geos_base::notice(const char* msg, void* data) {
    //std::cout << "geos notice " << msg << std::endl;
    //lastnotice=std::string(msg);
    
}

std::string geos_base::version() {
    return std::string(GEOSversion());
}

GEOSContextHandle_t geos_base::handle() { return handle_; }
        
//std::string geos_base::lastnotice="";
//std::string geos_base::lasterror="";
void* geos_base::data=nullptr;



coords convert_coords(std::shared_ptr<geos_base> base, const GEOSGeometry* geom) {
    const GEOSCoordSequence* cs = GEOSGeom_getCoordSeq_r(base->handle(), geom);
    if (!cs) { throw std::domain_error("no coords"); }
    
    unsigned int sz=0;
    if (!GEOSCoordSeq_getSize_r(base->handle(), cs, &sz)) {
        throw std::domain_error("size failed");
    }
    
    coords result(sz);
    if (sz==0) { return result; }
    
    for (size_t i=0; i < sz; i++) {
        if (!GEOSCoordSeq_getX_r(base->handle(), cs, i, &result[i].first)) { throw std::domain_error("size failed"); }
        if (!GEOSCoordSeq_getY_r(base->handle(), cs, i, &result[i].second)) { throw std::domain_error("size failed"); }
    }
    
    
    return result;
}

std::vector<coords> convert_rings(std::shared_ptr<geos_base> base, GEOSGeometry* geom) {
    
    
    int nr = GEOSGetNumInteriorRings_r(base->handle(), geom);
    if (nr<0) { throw std::domain_error("not a polygon? [type="+std::string(GEOSGeomType_r(base->handle(), geom))+"] (GetNumInteriorRings)"); }
    
    std::vector<coords> result(1+nr);
    
    const GEOSGeometry* ring;
    ring = GEOSGetExteriorRing_r(base->handle(), geom);
    if (!ring) { throw std::domain_error("not a polygon? [type="+std::string(GEOSGeomType_r(base->handle(), geom))+"] (GetExteriorRing)"); }
    
    result[0] = convert_coords(base, ring);
    
    for (int i=0; i < nr; i++) {
        ring = GEOSGetInteriorRingN_r(base->handle(), geom, i);
        if (!ring) { throw std::domain_error("not a polygon? [type="+std::string(GEOSGeomType_r(base->handle(), geom))+"] (GetInteriorRing "+std::to_string(i)+"/"+std::to_string(nr)+")"); }
        result[i+1] = convert_coords(base, ring);
    }
    return result;
}

geos_geometry::geos_geometry(std::shared_ptr<geos_base> base_, GEOSGeometry* geom_) : base(base_), geom(geom_), bounds_set(false), bounds_{0,0,0,0} {}

int geos_geometry::type_id() {
    return GEOSGeomTypeId_r(base->handle(), geom);
}

std::string geos_geometry::type_() {
    return std::string(GEOSGeomType_r(base->handle(), geom));
}

std::string geos_geometry::wkt() {
    char* c = GEOSGeomToWKT_r(base->handle(), geom);
    std::string s(c);
    free(c);
    return s;
}
double geos_geometry::area() {
    double a;
    GEOSArea_r(base->handle(), geom, &a);
    return a;
}

double geos_geometry::length() {
    double a;
    GEOSLength_r(base->handle(), geom, &a);
    return a;
}

coords geos_geometry::get_coords() {
    return convert_coords(base, geom);
}

std::vector<coords> geos_geometry::get_rings() {
    return convert_rings(base, geom);
}

std::vector<std::shared_ptr<geos_geometry>> geos_geometry::get_geometries() {
    
    int ng = GEOSGetNumGeometries_r(base->handle(), geom);
    std::vector<std::shared_ptr<geos_geometry>> result;
    for (int i=0; i < ng; i++) {
        const GEOSGeometry* g = GEOSGetGeometryN_r(base->handle(), geom, i);
        if (!g) { continue; }
        GEOSGeometry* gc = GEOSGeom_clone_r(base->handle(), g);
        if (!gc) { continue; }
        result.push_back(std::make_shared<geos_geometry>(base, gc));
    }
    return result;
}

int geos_geometry::get_geometries_raw(std::vector<GEOSGeometry*>& collection) {
    int ng = GEOSGetNumGeometries_r(base->handle(), geom);
    
    for (int i=0; i < ng; i++) {
        const GEOSGeometry* g = GEOSGetGeometryN_r(base->handle(), geom, i);
        if (!g) { continue; }
        GEOSGeometry* gc = GEOSGeom_clone_r(base->handle(), g);
        if (!gc) { continue; }
        collection.push_back(gc);
    }
    return ng;
}

std::shared_ptr<geos_geometry> geos_geometry::simplify(double tol) {
    GEOSGeometry* result = GEOSTopologyPreserveSimplify_r(base->handle(), geom, tol);
    if (!result) { throw std::domain_error("failed"); }
    return std::make_shared<geos_geometry>(base, result);
}

std::shared_ptr<geos_geometry> geos_geometry::buffer(double width, int quadsegs, int endCapStyle, int joinStyle, double mitreLimit) {
    GEOSGeometry* result = GEOSBufferWithStyle_r(base->handle(), geom, 
        width, quadsegs, endCapStyle, joinStyle, mitreLimit);
    
    if (!result) { throw std::domain_error("failed"); }
    return std::make_shared<geos_geometry>(base, result);
}

std::shared_ptr<geos_geometry> geos_geometry::buffer_zero() {
    return buffer(0, 16, 1, 1, 0);
}

bool geos_geometry::intersects(std::shared_ptr<geos_geometry> other) {
    char r = GEOSIntersects_r(base->handle(), geom, other->get_geometry());
    return r==1;
}

bool geos_geometry::equals_exact(std::shared_ptr<geos_geometry> other, double tol) {
    char r = GEOSEqualsExact_r(base->handle(), geom, other->get_geometry(), tol);
    return r==1;
}

bool geos_geometry::equals(std::shared_ptr<geos_geometry> other) {
    char r = GEOSEquals_r(base->handle(), geom, other->get_geometry());
    return r==1;
}

bool geos_geometry::disjoint(std::shared_ptr<geos_geometry> other) {
    char r = GEOSDisjoint_r(base->handle(), geom, other->get_geometry());
    return r==1;
}
bool geos_geometry::contains(std::shared_ptr<geos_geometry> other) {
    char r = GEOSContains_r(base->handle(), geom, other->get_geometry());
    return r==1;
}


std::shared_ptr<geos_geometry> geos_geometry::intersection(std::shared_ptr<geos_geometry> other) {
    GEOSGeometry* result = GEOSIntersection_r(base->handle(), geom, other->get_geometry());
    if (!result) { throw std::domain_error("failed " + base->lasterror); }
    return std::make_shared<geos_geometry>(base, result);
}

std::shared_ptr<geos_geometry> geos_geometry::clip(double a, double b, double c, double d) {
    if ((a>c) || (b>d)) {
        throw std::domain_error("clip box invalid");
    }
    GEOSGeometry* result = GEOSClipByRect_r(base->handle(), geom, a,b,c,d);
    if (!result) { throw std::domain_error("failed " + base->lasterror); }
    return std::make_shared<geos_geometry>(base, result);
}


std::shared_ptr<geos_geometry> geos_geometry::unary_union() {
    GEOSGeometry* result = GEOSUnaryUnion_r(base->handle(), geom);
    if (!result) { throw std::domain_error("failed " + base->lasterror); }
    return std::make_shared<geos_geometry>(base, result);
}

std::shared_ptr<geos_geometry> geos_geometry::difference(std::shared_ptr<geos_geometry> other) {
    GEOSGeometry* result = GEOSDifference_r(base->handle(), geom, other->get_geometry());
    if (!result) { throw std::domain_error("failed " + base->lasterror); }
    return std::make_shared<geos_geometry>(base, result);
}
std::shared_ptr<geos_geometry> geos_geometry::union_(std::shared_ptr<geos_geometry> other) {
    GEOSGeometry* result = GEOSUnion_r(base->handle(), geom, other->get_geometry());
    if (!result) { throw std::domain_error("failed " + base->lasterror); }
    return std::make_shared<geos_geometry>(base, result);
}

std::shared_ptr<geos_geometry> geos_geometry::envelope() {
    GEOSGeometry* result = GEOSEnvelope_r(base->handle(), geom);
    if (!result) { throw std::domain_error("failed " + base->lasterror); }
    return std::make_shared<geos_geometry>(base, result);
}

bool geos_geometry::is_valid() {
    return GEOSisValid_r(base->handle(), geom);
}

std::string geos_geometry::is_valid_reason() {
    return std::string(GEOSisValidReason_r(base->handle(), geom));
}

bool geos_geometry::is_empty() {
    return GEOSisEmpty_r(base->handle(), geom);
}

std::shared_ptr<geos_geometry> geos_geometry::linemerge() {
    GEOSGeometry* result = GEOSLineMerge_r(base->handle(), geom);
    if (!result) { throw std::domain_error("failed " + base->lasterror); }
    return std::make_shared<geos_geometry>(base, result);
}

std::shared_ptr<geos_geometry> geos_geometry::boundary() {
    GEOSGeometry* result = GEOSBoundary_r(base->handle(), geom);
    if (!result) { throw std::domain_error("failed " + base->lasterror); }
    return std::make_shared<geos_geometry>(base, result);
}

std::shared_ptr<geos_geometry> geos_geometry::centroid() {
    GEOSGeometry* result = GEOSGetCentroid_r(base->handle(), geom);
    if (!result) { throw std::domain_error("failed " + base->lasterror); }
    return std::make_shared<geos_geometry>(base, result);
}

    

std::string geos_geometry::wkb() {
    size_t sz;
    unsigned char* c = GEOSGeomToWKB_buf_r(base->handle(), geom, &sz);
    std::string s(reinterpret_cast<const char*>(c), sz);
    free(c);
    return s;
}

std::tuple<double,double,double,double> geos_geometry::bounds() {
    if (!bounds_set) {
        calc_bounds();
        bounds_set=true;
    }
    return bounds_;
}


void geos_geometry::calc_bounds() {
    if (is_empty()) {
        return;
    }
    if (type_id()==0) {
        auto p = get_coords();
        if (p.size()!=1) { throw std::domain_error("???"); }
        auto bl = p[0];
        bounds_=std::make_tuple(bl.first,bl.second,bl.first,bl.second);
        return;
    }

    auto env = envelope();
    if (env->is_empty()) {
        return;
    }
    
    if (env->type_id()==0) {
        auto p = env->get_coords();
        if (p.size()!=1) { throw std::domain_error("???"); }
        auto bl = p[0];
        bounds_=std::make_tuple(bl.first,bl.second,bl.first,bl.second);
        return;
    }
        
    try {
        
        auto cc = env->get_rings();
        if (cc.size()!=1) { throw std::domain_error("???"); }
        if (cc[0].size()!=5) { throw std::domain_error("???"); }
    
        auto bl = cc[0][0];
        auto tr = cc[0][2];
        bounds_=std::make_tuple(bl.first,bl.second,tr.first,tr.second);
    } catch (std::exception& ex) {
        throw std::domain_error("calc_bounds ["s+type_()+"]: "s+ex.what());
    }
}

geos_geometry::~geos_geometry() {
    GEOSGeom_destroy(geom);
}

GEOSGeometry* geos_geometry::get_geometry() { return geom; }
std::shared_ptr<geos_base> geos_geometry::get_base() { return base; }
    

py::object get_coords_py(std::shared_ptr<geos_geometry> geom) {
    if (geom->is_empty()) {
        return py::none();
    }
    
    int tid=geom->type_id();
    if (geom->type_id() == 0) {
        return py::make_tuple(geom->type_(), py::cast(geom->get_coords()[0]));
    }
    if (geom->type_id() == 1) {
        return py::make_tuple(geom->type_(), py::cast(geom->get_coords()));
    }
    if (geom->type_id() == 2) {
        return py::make_tuple(geom->type_(), py::cast(geom->get_coords()));
    }
    
    if (geom->type_id() == 3) {
        return py::make_tuple(geom->type_(), py::cast(geom->get_rings()));
    }
    
    if ((tid==4) || (tid==5) || (tid==6) || (tid==7)) {
        py::list result;
        for (const auto& g: geom->get_geometries()) {
            result.append(get_coords_py(g));
        }
        return py::make_tuple(geom->type_(), result);
    }
    throw std::domain_error("unexpected type "+geom->type_());
    return py::none();
}

picojson::value geos_point_geojson(const std::pair<double,double>& pt, int dp) {
    return picojson::value(
        picojson::array{
                picojson::value::with_precision(pt.first,dp),
                picojson::value::with_precision(pt.second,dp)
        }
    );
}

picojson::value geos_ring_geojson(const coords& ccs, int dp) {
    picojson::array rr; rr.reserve(ccs.size());
    for (const auto& c: ccs) {
        rr.push_back(geos_point_geojson(c, dp));
    }
    return picojson::value(rr);
}

picojson::value geos_geometry_geojson(std::shared_ptr<geos_geometry> geom, int dp) {
    picojson::object result;
    result["type"] = picojson::value(geom->type_());
    
    
    if (geom->type_id()==7) {
        picojson::array gg;
        for (auto g: geom->get_geometries()) {
            gg.push_back(geos_geometry_geojson(g,dp));
        }
        result["geometries"] = picojson::value(gg);
        return picojson::value(result);
    }
    
    
    
    if (geom->type_id()==0) {
        auto cc=geom->get_coords();
        if (cc.size()!=1) { throw std::domain_error("??"); }
        result["coordinates"] = geos_point_geojson(cc[0],dp);
    } else if (geom->type_id()==1) {
        auto cc=geom->get_coords();
        result["coordinates"] = geos_ring_geojson(cc,dp);
    } else if (geom->type_id()==2) {
        throw std::domain_error("geojson of linearring ??");
    } else if (geom->type_id()==3) {
        auto rings=geom->get_rings();
        picojson::array pp;
        for (auto rr: rings) {
            pp.push_back(picojson::value(geos_ring_geojson(rr,dp)));
        }
        result["coordinates"]=picojson::value(pp);
    } else if (geom->type_id()==4) {
        
        picojson::array mp;
        for (auto gg: geom->get_geometries()) {
            if (gg->type_id()!=0) { throw std::domain_error("multipoint contains a "+gg->type_()); }
            
            auto cc=gg->get_coords();
            if (cc.size()!=1) { throw std::domain_error("??"); }
            mp.push_back(picojson::value(geos_point_geojson(cc[0],dp)));
        }
        result["coordinates"]=picojson::value(mp);
    } else if (geom->type_id()==5) {
        
        picojson::array ml;
        for (auto gg: geom->get_geometries()) {
            if (gg->type_id()!=1) { throw std::domain_error("multiline contains a "+gg->type_()); }
            
            auto cc=gg->get_coords();
            
            ml.push_back(picojson::value(geos_ring_geojson(cc,dp)));
        }
        result["coordinates"]=picojson::value(ml);
    } else if (geom->type_id()==6) {
        
        picojson::array mp;
        for (auto gg: geom->get_geometries()) {
            if (gg->type_id()!=3) { throw std::domain_error("multipolygon contains a "+gg->type_()); }
            
            auto rings=gg->get_rings();
            picojson::array pp;
            for (auto rr: rings) {
                pp.push_back(picojson::value(geos_ring_geojson(rr,dp)));
            }
            
            mp.push_back(picojson::value(pp));
        }
        result["coordinates"]=picojson::value(mp);
    }
    
    return picojson::value(result);
}

std::shared_ptr<geos_geometry> geos_geometry_create_wkb(std::shared_ptr<geos_base> base, const std::string& wkb) {
    GEOSGeometry* geom = GEOSGeomFromWKB_buf_r(base->handle(), reinterpret_cast<const unsigned char*>(wkb.c_str()), wkb.size());
    if (!geom) { throw std::domain_error("no geom created"); }
    return std::make_shared<geos_geometry>(base, geom);
}

std::shared_ptr<geos_geometry> geos_geometry_create_point(std::shared_ptr<geos_base> base, double x, double y) {
    GEOSCoordSequence* cs = GEOSCoordSeq_create_r(base->handle(), 1, 2);
    GEOSCoordSeq_setX_r(base->handle(), cs, 0, x);
    GEOSCoordSeq_setY_r(base->handle(), cs, 0, y);
    GEOSGeometry* geom = GEOSGeom_createPoint_r(base->handle(), cs);
    if (!geom) { throw std::domain_error("no geom created"); }
    return std::make_shared<geos_geometry>(base, geom);
}

GEOSCoordSequence* make_coordseq(const std::shared_ptr<geos_base>& base, const coords& ccs) {
    GEOSCoordSequence* cs = GEOSCoordSeq_create_r(base->handle(), ccs.size(), 2);
    for (size_t i=0; i < ccs.size(); i++) {
        GEOSCoordSeq_setX_r(base->handle(), cs, i, ccs[i].first);
        GEOSCoordSeq_setY_r(base->handle(), cs, i, ccs[i].second);
    }
    return cs;
}

std::shared_ptr<geos_geometry> geos_geometry_create_linestring(std::shared_ptr<geos_base> base, const coords& ccs) {
    GEOSCoordSequence* cs = make_coordseq(base, ccs);
    GEOSGeometry* geom = GEOSGeom_createLineString_r(base->handle(), cs);
    if (!geom) { throw std::domain_error("no geom created"); }
    return std::make_shared<geos_geometry>(base, geom);
}

std::shared_ptr<geos_geometry> geos_geometry_create_polygon(std::shared_ptr<geos_base> base, const std::vector<coords>& rings_in) {
    if (rings_in.empty()) {
        return std::make_shared<geos_geometry>(base, GEOSGeom_createEmptyPolygon_r(base->handle()));
    }
    
    std::vector<GEOSGeometry*> rings(rings_in.size());
    for (size_t i=0; i < rings_in.size(); i++) {
        rings[i] = GEOSGeom_createLinearRing_r(base->handle(), make_coordseq(base, rings_in[i]));
    }
    
    GEOSGeometry* geom = GEOSGeom_createPolygon_r(base->handle(), rings[0], &rings[1], rings.size()-1);
    if (!geom) { throw std::domain_error("no geom created"); }
    return std::make_shared<geos_geometry>(base, geom);
}
int as_mm_typeid(int ti) {
    if (ti==0) { return 4; }
    if (ti==1) { return 5; }
    if (ti==2) { throw std::domain_error("multi of rings??"); }
    if (ti==3) { return 6; }
    return ti;
}
std::shared_ptr<geos_geometry> geos_geometry_create_collection(std::shared_ptr<geos_base> base, const std::vector<std::shared_ptr<geos_geometry>>& geoms_in) {
    if (geoms_in.empty()) {
        return std::make_shared<geos_geometry>(base, GEOSGeom_createEmptyCollection_r(base->handle(), 7));
    }
    
    int type_=-1;
    std::vector<GEOSGeometry*> geoms;
    for (auto& g: geoms_in) {
        int gt=as_mm_typeid(g->type_id());
        if (type_==-1) { type_=gt; }
        else if (type_!=gt) { type_=7; }
        
        g->get_geometries_raw(geoms);
    }
    if (geoms.empty()) {
        return std::make_shared<geos_geometry>(base, GEOSGeom_createEmptyCollection_r(base->handle(),7));
    }
    if (geoms.size()==1) {
        return std::make_shared<geos_geometry>(base, geoms[0]);
    }
    
    GEOSGeometry* result = GEOSGeom_createCollection_r(base->handle(), type_, &geoms[0], geoms.size());
    if (!result) { throw std::domain_error("no geom created"); }
    return std::make_shared<geos_geometry>(base, result);
}

std::shared_ptr<geos_geometry> filter_collection(int type_id, std::shared_ptr<geos_geometry> src) {
    auto base=src->get_base();
    if (src->is_empty()) {
        return src;
    }
    if (type_id==7) { throw std::domain_error("????"); }
    if (type_id==0) { type_id=4; }
    if (type_id==1) { type_id=5; }
    if (type_id==3) { type_id=6; }
           
        
    
    std::vector<GEOSGeometry*> geoms;
    
    
    for (auto& g: src->get_geometries()) {
        if (as_mm_typeid(g->type_id())==type_id) {
            g->get_geometries_raw(geoms);
        }
        
    }
    if (geoms.empty()) {
        return std::make_shared<geos_geometry>(base, GEOSGeom_createEmptyCollection_r(base->handle(),7));
    }
    if (geoms.size()==1) {
        return std::make_shared<geos_geometry>(base, geoms[0]);
    }
    
    GEOSGeometry* result = GEOSGeom_createCollection_r(base->handle(), type_id, &geoms[0], geoms.size());
    if (!result) { throw std::domain_error("no geom created"); }
    return std::make_shared<geos_geometry>(base, result);
}


std::shared_ptr<geos_geometry> geos_geometry_transform_point(std::shared_ptr<geos_geometry> src, const std::function<std::pair<double,double>(double,double)>& func) {
    auto cc = src->get_coords();
    if (cc.size()!=1) { throw std::domain_error("???"); }
    
    auto p = func(cc[0].first, cc[0].second);
    return geos_geometry_create_point(src->get_base(), p.first, p.second);
}

void geos_geometry_transform_line(coords& cc, const std::function<std::pair<double,double>(double,double)>& func) {
    if (cc.empty()) { return; }
    for( size_t i =0; i < cc.size(); i++) {
        cc[i] = func(cc[i].first, cc[i].second);
    }
}

std::shared_ptr<geos_geometry> geos_geometry_transform_linestring(std::shared_ptr<geos_geometry> src, const std::function<std::pair<double,double>(double,double)>& func) {
    auto cc = src->get_coords();
    
    geos_geometry_transform_line(cc, func);
    
    return geos_geometry_create_linestring(src->get_base(), cc);
}

std::shared_ptr<geos_geometry> geos_geometry_transform_polygon(std::shared_ptr<geos_geometry> src, const std::function<std::pair<double,double>(double,double)>& func) {
    auto rr = src->get_rings();
    
    for (auto& cc: rr) {
        geos_geometry_transform_line(cc, func);
    }
    return geos_geometry_create_polygon(src->get_base(), rr);
}




std::shared_ptr<geos_geometry> geos_geometry_transform(std::shared_ptr<geos_geometry> src, const std::function<std::pair<double,double>(double,double)>& func) {
    
    if (src->type_id() == 0) {
        return geos_geometry_transform_point(src, func);
    } else if (src->type_id() == 1) {
        return geos_geometry_transform_linestring(src, func);
    } else if (src->type_id() == 2) {
        throw std::domain_error("?? linearring");
        
    } else if (src->type_id() == 3) {
        return geos_geometry_transform_polygon(src, func);
    }
    
    auto src_geoms = src->get_geometries();
    std::vector<std::shared_ptr<geos_geometry>> out_geoms;
    out_geoms.reserve(src_geoms.size());
    for (auto g: src_geoms) {
        out_geoms.push_back(geos_geometry_transform(g, func));
    }
    return geos_geometry_create_collection(src->get_base(), out_geoms);
}
    

std::shared_ptr<geos_geometry> geos_geometry_transform_scale(std::shared_ptr<geos_geometry> src, double sc) {
    return geos_geometry_transform(src, [sc](double x, double y) { return std::make_pair(x*sc, y*sc); });
}

std::shared_ptr<geos_geometry> geos_geometry_transform_translate(std::shared_ptr<geos_geometry> src, double x0, double y0) {
    return geos_geometry_transform(src, [x0,y0](double x, double y) { return std::make_pair(x+x0, y+y0); });
}

std::shared_ptr<geos_geometry> geos_geometry_transform_translate_scale(std::shared_ptr<geos_geometry> src, double x0, double y0, double sc) {
    return geos_geometry_transform(src, [sc, x0,y0](double x, double y) { return std::make_pair( (x+x0)*sc, (y+y0)*sc); });
}

std::shared_ptr<geos_geometry> geos_geometry_transform_tile(std::shared_ptr<geos_geometry> src, int64 tx, int64 ty, int64 tz) {
    
    bbox bx = tile_bound(tx,ty,tz);
    double sc = (std::get<2>(bx)-std::get<0>(bx))/256.0;
    double x0 = std::get<0>(bx);
    double y0 = std::get<1>(bx);
    
    return geos_geometry_transform(src, [x0,y0,sc](double x, double y) { return std::make_pair( x0 + x*sc, y0+y*sc); });
}
    
    
    
    

    
std::shared_ptr<geos_geometry> geos_geometry_create_box(std::shared_ptr<geos_base> base, double a, double b, double c, double d) {
    
    coords cc{{a,b},{c,b},{c,d},{a,d},{a,b}};
    return geos_geometry_create_polygon(base, {cc});
}
  
void iter_geometry_tiles(std::shared_ptr<geos_geometry> geom,
    int64 x, int64 y, int64 z, int64 minzoom, int64 maxzoom,
    const std::function<void(int64,int64,int64)>& cb) {


    double a,b,c,d;
    std::tie(a,b,c,d) = tile_bound(x,y,z,0);
    
    auto gc = geom->clip(a,b,c,d);
    
    if (!gc) {
        std::cout << "clip failed " << x << " " << y << " " << z << std::endl;
        throw std::domain_error("??");
    }
    if (gc->is_empty()) {
        return;
    }
    
    if (z>=minzoom) {
        cb(x,y,z);
    }
    if (z >= maxzoom) { return; }
    
    iter_geometry_tiles(gc, 2*x,  2*y,  z+1,minzoom, maxzoom,cb);
    iter_geometry_tiles(gc, 2*x+1,2*y,  z+1,minzoom, maxzoom,cb);
    iter_geometry_tiles(gc, 2*x,  2*y+1,z+1,minzoom, maxzoom,cb);
    iter_geometry_tiles(gc, 2*x+1,2*y+1,z+1,minzoom, maxzoom,cb);
    
}



    
    
std::vector<std::tuple<int64,int64,int64>> geometry_tiles(
    std::shared_ptr<geos_geometry> geom, int64 x, int64 y, int64 z, int64 minzoom, int64 maxzoom) {
    
    std::vector<std::tuple<int64,int64,int64>> result;
    std::function<void(int64,int64,int64)> cb =
        [&result](int64 xx, int64 yy, int64 zz) { result.push_back(std::make_tuple(xx,yy,zz)); };
    
    
    iter_geometry_tiles(geom, x, y, z, minzoom, maxzoom, cb);
    return result;
};

xyz next_branch(const xyz& t) {
    int64 x,y,z;
    std::tie(x,y,z)=t;
    
    if (z==0) {
        return std::make_tuple(-1,-1,-1);
    }
    
    if ((x%2)==0) {
        return std::make_tuple(x+1,y,z);
    }
    if ((y%2)==0) {
        return std::make_tuple(x-1,y+1,z);
    }
    
    return next_branch(std::make_tuple(x/2,y/2,z-1));
}

xyz next_tile(const xyz& t, int64 maxzoom) {
    
    if (std::get<2>(t) < maxzoom) {
        return std::make_tuple(std::get<0>(t)*2, std::get<1>(t)*2, std::get<2>(t)+1);
    }
    return next_branch(t);
}
    

GeometryTilesIter::GeometryTilesIter(std::shared_ptr<geos_geometry> geom_,
    std::tuple<int64,int64,int64> min_tile_,
    std::tuple<int64,int64,int64> max_tile_,
    int64 maxzoom_) : geom(geom_), min_tile(min_tile_), max_tile(max_tile_), maxzoom(maxzoom_), curr(min_tile_) {}
    

bool check_poly_tile(const std::shared_ptr<geos_geometry>& geom, const xyz& t) {
    if (!geom) { return true; }
    bbox bx = tile_bound(std::get<0>(t),std::get<1>(t),std::get<2>(t),0);
    try {
        auto gc = geom->clip(std::get<0>(bx),std::get<1>(bx),std::get<2>(bx),std::get<3>(bx));
        return !gc->is_empty();
    } catch(std::exception& ex) {
        return false;
    }
    
    return false;
}

xyz tile_round(const xyz& t, int64 zoom) {
    if (std::get<2>(t) <= zoom) { return t; }
    
    int64 zz = 1 << (std::get<2>(t)-zoom);
    return xyz{std::get<0>(t)/zz, std::get<1>(t)/zz, zoom};
}


int64 tile_square(const xyz& t, int64 l) {
    if (l >= std::get<2>(t)) { throw std::domain_error("out of range"); }
    
    int64 xx,yy,zz;
    std::tie(xx,yy,zz) = t;
    zz = (zz-l-1);
    xx = (xx>>zz)%2;
    yy = (yy>>zz)%2;
    
    return xx + 2*yy;
}
    
std::string tile_string(const xyz& t) {
    size_t z = std::get<2>(t);
    std::string r('\0', z);
    for (size_t i=0; i < z; i++) {
        r[i] = 'A'+tile_square(t, i);
    }
    return r;
}
std::ostream& operator<<(std::ostream& os, const xyz& t) {
    os << "[" << std::get<0>(t) << ", " << std::get<1>(t) << ", " << std::get<2>(t) << "]";
    return os;
}
bool tile_lessthan(const xyz& left, const xyz& right) {
    if (left==right) { return false; }
    if (std::get<2>(right)<0) { return true; }
    if (std::get<2>(left)<0) {
        std::cout << "tile_lessthan " << left << ", " << right << std::endl;        
        throw std::domain_error("wtf");
    }
        
    
    size_t mz = std::min(std::get<2>(left), std::get<2>(right));
    
    if (mz==0) { return std::get<2>(left) < std::get<2>(right); }
    
    for (size_t i=0; i < mz; i++) {
        int64 l = tile_square(left,i);
        int64 r = tile_square(right,i);
        if (l<r) { return true; }
        if (l>r) { return false; }
    }
    
    return std::get<2>(left)<std::get<2>(right);
}

std::tuple<int64,int64,int64> GeometryTilesIter::next() {
    
    while ((std::get<0>(curr)>=0) && tile_lessthan(curr, max_tile)) {
        if (check_poly_tile(geom, curr)) {
            xyz ans=curr;
            curr = next_tile(curr, maxzoom);
            return ans;
        } else {
            curr=next_branch(curr);
            
        }
    }
    
    return std::make_tuple(-1,-1,-1);
}

std::vector<xyz> iter_all(GeometryTilesIter& gti) {
    std::vector<xyz> result;
    auto t = gti.next();
    while (std::get<2>(t)>=0) {
        result.push_back(t);
        t = gti.next();
    }
    return result;
}   
    
std::vector<xyz> geometry_tiles_all_between(std::shared_ptr<geos_geometry> geom, const xyz& mintile, const xyz& maxtile, int64 maxzoom) {
    GeometryTilesIter gti(geom,mintile,maxtile,maxzoom);
    
    return iter_all(gti);
}   
    
void export_geos_geometry(py::module& m) {
    py::class_<geos_base, std::shared_ptr<geos_base>>(m,"geos_base")
        .def(py::init<>())
        .def("version", &geos_base::version)
        .def("wkb", &geos_geometry_create_wkb)
        .def("point", &geos_geometry_create_point)
        .def("linestring", &geos_geometry_create_linestring)
        .def("polygon", &geos_geometry_create_polygon)
        .def("collection", &geos_geometry_create_collection)
        .def("box", &geos_geometry_create_box)
    ;
    py::class_<geos_geometry, std::shared_ptr<geos_geometry>>(m, "geos_geometry")
        .def_property_readonly("type", &geos_geometry::type_)
        .def_property_readonly("type_id", &geos_geometry::type_id)
        .def_property_readonly("area", &geos_geometry::area)
        .def_property_readonly("length", &geos_geometry::length)
        .def_property_readonly("wkt", &geos_geometry::wkt)
        .def_property_readonly("wkb", [](std::shared_ptr<geos_geometry> g) { return py::bytes(g->wkb()); })
        
        .def_property_readonly("mapping", get_coords_py)
        .def_property_readonly("is_valid", &geos_geometry::is_valid)
        .def_property_readonly("is_valid_reason", &geos_geometry::is_valid_reason)
        .def_property_readonly("is_empty", &geos_geometry::is_empty)
        
        .def("buffer", &geos_geometry::buffer, py::arg("buffer"), py::arg("quadsegs")=16, py::arg("endcap")=1, py::arg("join")=1, py::arg("mitre")=0.0)
        .def("simplify", &geos_geometry::simplify)
        
        .def("intersects", &geos_geometry::intersects)
        .def("equals", &geos_geometry::equals)
        .def("equals_exact", &geos_geometry::equals_exact)
        .def("disjoint", &geos_geometry::disjoint)
        .def("contains", &geos_geometry::contains)
        
        .def("difference", &geos_geometry::difference)
        .def("union", &geos_geometry::union_)
        .def("unary_union", &geos_geometry::unary_union)
        .def("intersection", &geos_geometry::intersection)
        .def("clip", &geos_geometry::clip)
        .def("envelope", &geos_geometry::envelope)
        
        .def("boundary", &geos_geometry::boundary)
        .def("linemerge", &geos_geometry::linemerge)
        .def("centroid", &geos_geometry::centroid)
        
        .def_property_readonly("bounds", &geos_geometry::bounds)
        
        .def("geojson", [](std::shared_ptr<geos_geometry> geom, int dp) { return py::bytes(geos_geometry_geojson(geom, dp).serialize()); })
        .def("get_geometries", &geos_geometry::get_geometries)
        
    ;   
    
    m.def("geometry_tiles", &geometry_tiles);
    
    m.def("transform", &geos_geometry_transform);
    m.def("transform_scale", &geos_geometry_transform_scale);
    m.def("transform_translate", &geos_geometry_transform_translate);
    m.def("transform_translate_scale", &geos_geometry_transform_translate_scale);
    m.def("transform_tile", &geos_geometry_transform_tile);
    
    
    py::class_<GeometryTilesIter>(m, "GeometryTilesIter")
        .def(py::init<std::shared_ptr<geos_geometry>,xyz,xyz,int64>())
        .def("next", &GeometryTilesIter::next)
        
        .def("iter_all", &iter_all)
    ;
    
    m.def("tile_round", &tile_round);
    m.def("next_tile", &next_tile);
    m.def("next_branch", &next_branch);
    m.def("tile_square", &tile_square);
    m.def("tile_string", &tile_string);
    m.def("tile_lessthan", &tile_lessthan);
}


