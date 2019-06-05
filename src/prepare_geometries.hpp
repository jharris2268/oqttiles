#ifndef PREPARE_GEOMETRIES_HPP
#define PREPARE_GEOMETRIES_HPP

#include "geos_wrapper.hpp"
#include <oqt/elements/block.hpp>
#include <memory>

typedef std::tuple<int64,int64,int64> xyz;
typedef std::map<xyz,std::string> packed_tiles;
typedef std::map<std::string, picojson::value> property_map;
typedef std::tuple<std::string,std::string,property_map,bool> feature_table_entry;
typedef std::map<int64,std::vector<feature_table_entry>> feature_tables;


struct alt_feature_data {
    int64 id;
    std::string properties;
    double areaorlen;
    size_t geom_type;
    int64 np;
    std::string geom_data;
    int64 minzoom;
};

struct alt_tile_data {
    alt_tile_data() : tile_xyz(0,0,0) {}
    xyz tile_xyz;
    std::map<std::pair<xyz, std::string>, std::vector<alt_feature_data>> objs;
};


typedef std::map<
    std::tuple<int,std::string,std::string>,
    std::pair<int64,std::vector<std::string>>
    > feature_spec;

typedef std::pair<
    int64,        //geometry type
    std::vector<
        std::tuple<
            std::vector<std::pair<std::string,std::string>>, //rules
            int64,                                           //minzoom
            std::string                                      //tag
        >
    >> extra_tags_spec;




struct alt_tile_data_vec {
    alt_tile_data_vec() {}
    alt_tile_data_vec(xyz tile_xyz_, std::vector<std::shared_ptr<alt_tile_data>> tiles_, std::vector<xyz> to_finish_)
            : tile_xyz(tile_xyz_), tiles(tiles_), to_finish(to_finish_) {}
    xyz tile_xyz;
    std::vector<std::shared_ptr<alt_tile_data>> tiles;
    std::vector<xyz> to_finish;
};

typedef std::function<void(std::shared_ptr<alt_tile_data_vec>)> alt_tile_data_vec_cb;
std::function<void(std::shared_ptr<alt_tile_data>)> make_addblockstree_alt_cb(std::shared_ptr<geos_geometry> geom, alt_tile_data_vec_cb cb, int64 maxtile);

alt_tile_data_vec_cb make_addotherfeatures(std::shared_ptr<alt_tile_data> otherfeatures, alt_tile_data_vec_cb cb);

//std::function<void(std::shared_ptr<alt_tile_data_vec>)> make_collectaltblocks_cb(std::shared_ptr<geos_geometry> geom, std::vector<alt_tile_data_vec_cb> cb, int64 maxdepth);
#endif
