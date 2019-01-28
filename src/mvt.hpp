#ifndef MVT_HPP
#define MVT_HPP

#include "prepare_geometries.hpp"



alt_feature_data prep_alt_feature_data(int64 id, const property_map& pm, std::shared_ptr<geos_geometry> geom, int64 minzoom, int64 np);

std::shared_ptr<geos_geometry> read_geometry(std::shared_ptr<geos_base>, size_t geom_ty, const std::string& geom_data, int64 np);

alt_tile_data_vec_cb make_merge_interim_tiles_callback(bool compress, int64 maxlevel, bool sortobjs,  bool mergefeats, std::function<void(std::shared_ptr<packed_tiles>)> cb);
#endif
