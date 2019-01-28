#ifndef MAKETILEDATA_HPP
#define MAKETILEDATA_HPP

#include "prepare_geometries.hpp"



std::function<void(oqt::PrimitiveBlockPtr)> make_maketiledata_alt_callback(
    feature_spec features, std::map<std::string,extra_tags_spec> extra_tags,
    int64 max_tile, bool use_obj_tile, bbox filter_box, bool fix_geoms, std::function<void(std::shared_ptr<alt_tile_data>)> cb);

std::function<void(oqt::PrimitiveBlockPtr)> make_maketiledata_groupalt_callback(
feature_spec features, std::map<std::string,extra_tags_spec> extra_tags,
    int64 max_tile, bool use_obj_tile,
    bbox filter_box, bool fix_geoms, const std::vector<xyz>& qts, std::function<void(std::shared_ptr<alt_tile_data_vec>)> cb) ;
std::string tup_str(const xyz& t);

#endif
