#ifndef GEOS_WRAPPER_HPP
#define GEOS_WRAPPER_HPP

#include <picojson_prec.h>

#include "geos_c.h"

#include <string>
#include <vector>
#include <memory>
#include <atomic>
#include <functional>

typedef int64_t int64;
typedef std::tuple<double,double,double,double> bbox;
bool bbox_overlaps(const bbox& l, const bbox& r);

typedef std::tuple<int64,int64,int64> xyz;
std::ostream& operator<<(std::ostream& os, const xyz& t);

bbox tile_bound(int64 x, int64 y, int64 z, double buf=0);
double zoom_scale(int64 z);

class geos_base {
    public:
        geos_base();
        
        virtual ~geos_base();
        
        static void error(const char* msg, void* data);
        static void notice(const char* msg, void* data);
        
        std::string version();
        GEOSContextHandle_t handle();
        
        //static std::atomic<std::string> lastnotice;
        //static std::atomic<std::string> lasterror;
        std::string lasterror;
    private:
        GEOSContextHandle_t handle_;
        
        static void* data;
};
typedef std::vector<std::pair<double,double>> coords;


class geos_geometry {
    public:
        geos_geometry(std::shared_ptr<geos_base> base_, GEOSGeometry* geom_);
        
        int type_id();
        
        std::string type_();
        std::string wkt();
        double area();
        
        double length();
        
        coords get_coords();
        std::vector<coords> get_rings();
        
        std::vector<std::shared_ptr<geos_geometry>> get_geometries();
        
        int get_geometries_raw(std::vector<GEOSGeometry*>& collection);
        
        std::shared_ptr<geos_geometry> simplify(double tol);
        
        std::shared_ptr<geos_geometry> buffer(double width, int quadsegs, int endCapStyle, int joinStyle, double mitreLimit);
        
        std::shared_ptr<geos_geometry> buffer_zero();
        
        bool intersects(std::shared_ptr<geos_geometry> other);
        
        bool equals_exact(std::shared_ptr<geos_geometry> other, double tol);
        
        bool equals(std::shared_ptr<geos_geometry> other);
        
        bool disjoint(std::shared_ptr<geos_geometry> other);
        bool contains(std::shared_ptr<geos_geometry> other);
        
        
        std::shared_ptr<geos_geometry> intersection(std::shared_ptr<geos_geometry> other);
        std::shared_ptr<geos_geometry> clip(double a, double b, double c, double d);
        
        std::shared_ptr<geos_geometry> unary_union();
        
        std::shared_ptr<geos_geometry> difference(std::shared_ptr<geos_geometry> other);
        std::shared_ptr<geos_geometry> union_(std::shared_ptr<geos_geometry> other);
        std::shared_ptr<geos_geometry> envelope();
        
        bool is_valid();
        std::string is_valid_reason();
        bool is_empty();
        
        std::shared_ptr<geos_geometry> linemerge();
        
        std::shared_ptr<geos_geometry> boundary();
        
        std::shared_ptr<geos_geometry> centroid();
        
            
        
        std::string wkb();
        
        std::tuple<double,double,double,double> bounds();
        void calc_bounds();
        
        virtual ~geos_geometry();
        
        GEOSGeometry* get_geometry();
        std::shared_ptr<geos_base> get_base();
    
    private:
        std::shared_ptr<geos_base> base;
        GEOSGeometry* geom;
        
        bool bounds_set;
        std::tuple<double,double,double,double> bounds_;
};

picojson::value geos_geometry_geojson(std::shared_ptr<geos_geometry> geom, int dp);

std::shared_ptr<geos_geometry> geos_geometry_create_wkb(std::shared_ptr<geos_base> base, const std::string& wkb);

std::shared_ptr<geos_geometry> geos_geometry_create_point(std::shared_ptr<geos_base> base, double x, double y);
std::shared_ptr<geos_geometry> geos_geometry_create_linestring(std::shared_ptr<geos_base> base, const coords& ccs);
std::shared_ptr<geos_geometry> geos_geometry_create_polygon(std::shared_ptr<geos_base> base, const std::vector<coords>& rings_in);
std::shared_ptr<geos_geometry> geos_geometry_create_collection(std::shared_ptr<geos_base> base, const std::vector<std::shared_ptr<geos_geometry>>& geoms_in);
std::shared_ptr<geos_geometry> geos_geometry_create_box(std::shared_ptr<geos_base> base, double a, double b, double c, double d);

std::shared_ptr<geos_geometry> geos_geometry_transform(std::shared_ptr<geos_geometry> src, const std::function<std::pair<double,double>(double,double)>& func);
std::shared_ptr<geos_geometry> geos_geometry_transform_scale(std::shared_ptr<geos_geometry> src, double sc);
std::shared_ptr<geos_geometry> geos_geometry_transform_translate(std::shared_ptr<geos_geometry> src, double x0, double y0);
std::shared_ptr<geos_geometry> geos_geometry_transform_translate_scale(std::shared_ptr<geos_geometry> src, double x0, double y0, double sc);    


std::shared_ptr<geos_geometry> filter_collection(int type_id, std::shared_ptr<geos_geometry> src);

void iter_geometry_tiles(std::shared_ptr<geos_geometry> geom,
    int64 x, int64 y, int64 z, int64 minzoom, int64 maxzoom,
    const std::function<void(int64,int64,int64)>& cb);

class GeometryTilesIter {
    public:
        GeometryTilesIter(std::shared_ptr<geos_geometry> geom,
            std::tuple<int64,int64,int64> min_tile,
            std::tuple<int64,int64,int64> max_tile,
            int64 maxzoom);
        
        std::tuple<int64,int64,int64> next();
    
    private:
        std::shared_ptr<geos_geometry> geom;
        std::tuple<int64,int64,int64> min_tile;
        std::tuple<int64,int64,int64> max_tile;
        int64 maxzoom;
        
        std::tuple<int64,int64,int64> curr;
};
xyz tile_round(const xyz& t, int64 zoom);
xyz next_tile(const xyz& t, int64 maxzoom);
xyz next_branch(const xyz& t);
int64 tile_square(const xyz& t, int64 l);
std::string tile_string(const xyz& t);
bool tile_lessthan(const xyz& left, const xyz& right);

std::vector<xyz> geometry_tiles_all_between(std::shared_ptr<geos_geometry> geom, const xyz& mintile, const xyz& maxtile, int64 maxzoom); 
    

#endif
