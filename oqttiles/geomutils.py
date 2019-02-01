#-----------------------------------------------------------------------
#
# This file is part of oqtsqlite
#
# Copyright (C) 2018 James Harris
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
#-----------------------------------------------------------------------


from . import _oqttiles


gs = _oqttiles.geos_base()

collection = gs.collection
point = gs.point
linestring = gs.linestring
polygon = gs.polygon
box = gs.box
wkb = gs.wkb

#geometry_tiles = geos_geometry_tiles

def convert_geom(obj):
    
    if obj.Type==3:
        pt=obj.LonLat.transform
        return point(pt.x, pt.y)
    
    if obj.Type==4:
        pp = [(p.x,p.y) for p in (ll.transform for ll in obj.LonLats)]
        return linestring(pp)
        
    if obj.Type==5:
        pp = [(p.x,p.y) for p in (ll.transform for ll in obj.LonLats)]
        return polygon([pp])
    
    if obj.Type==6:
        ext = [(p.x,p.y) for p in (ll.transform for ll in obj.OuterLonLats)]
        iis = [[(p.x,p.y) for p in (ll.transform for ll in ii)] for ii in obj.InnerLonLats]
        
        return polygon([ext]+iis)
