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
