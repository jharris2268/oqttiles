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
from __future__ import print_function


from . import _oqttiles
from .geomutils import box, gs

import oqt
from oqtsqlite import mbtiles as mbt

import time, os

tuple_to_str = lambda t: oqt.elements.quadtree_string(oqt.elements.quadtree_from_tuple(*t))


def get_geoms(geomsfn, ll, i,pf,kw):
    gg=[]
    oqt.geometry.read_blocks_geometry(geomsfn, oqt.utils.addto(gg),ll, **kw)
    if not pf is None:
        for g in gg:
            g.FileProgress=i*pf
            i+=1

    return gg
    
def get_geoms_iter(geomsfn, ll,target, kw):
    pf=100.0/len(ll)
    if sum(c/1024./1024 for a,b,c in ll)<1.8*target:
        yield get_geoms(geomsfn, [b for a,b,c in ll], 0,pf,kw)
        return
    
    mm=[]
    tot=0
    for i,(a,b,c) in enumerate(ll):
        mm.append(b)
        tot += c/1024./1024.
        if tot>target:
            yield get_geoms(geomsfn, mm,i,pf,kw)
            mm=[]
            tot=0
    
    if mm:
        get_geoms(geomsfn, mm,i,pf,kw)
    

def iter_geoms(geomsfn, hh, poly,minzoom,target=None):
    
    target=target or 50
    ll=[(a,b,c) for a,b,c in hh if (minzoom is None or (a&31)<=minzoom) and (poly is None or poly(a))]
    print("%d locs" % len(ll))
    kw={}
    if not minzoom is None:
        kw['minzoom']=minzoom
    if not poly is None and not poly.bounds is None:
        
        kw['bbox']=poly.bounds
    
    
    
    return get_geoms_iter(geomsfn,ll,target,kw)


def prep_spec(polypoint, all_tags=False):
    extra_tags={}
    
    extra_tags['boundary'] = (1,[([],4,'boundary'),([],4,'admin_level'),([],4,'name'), ([],4,'name:en')])    
    if all_tags:
        extra_tags['point'] = (0,[([],0,'*')])
        extra_tags['line'] = (1,[([],0,'*')])
        extra_tags['highway'] = (1,[([],0,'*')])
        extra_tags['polygon'] = (2,[([],0,'*')])
        extra_tags['building'] = (2,[([],0,'*')])
        
    else:
        extra_tags['point'] = (0,[
            ([],4,'name'),
            ([],4,'name:en'),
            ([],4,'capital'),
            ([],4,'population'),
            ([],12,'ref'),
            ([],10,'access'),
            ([],10,'icao'),
            ([],10,'iata'),    
            ([],13,'*')])
        
        
        
        extra_tags['line'] = (1,[
            ([],10,'intermittent'),
            ([],10,'seasonal'),
            ([],12,'bridge'),
            ([],12,'tunnel'),
            ([],12,'name'),
            ([],12,'name:en'),
            ([],13,'*')])
            
        
        
        extra_tags['building'] = (2,[
            ([('amenity','place_of_worship')], 12, 'amenity'),
            ([('building','train_station')],12,'building'),
            ([('aerialway','station')],12, 'aerialway'),
            ([('public_transport','station')],12,'public_transport'),
            ([('aeroway','terminal')],12, 'aeroway'),
            ([],13,'*'),
        ])
        
        
        extra_tags['polygon'] = (2,[
            ([],5,'natural'),
            ([],5,'landuse'),
            ([],5,'wetland'),
            ([],5,'surface'),
            ([('boundary','national_park')],5,'name'),
            ([('landuse','*')],5,'name'),
            ([('leisure','nature_reserve')],5,'name'),
            ([('military','danger_area')],5,'name'),
            ([('natural','*')],5,'name'),
            ([('place','island')],5,'name'),
            ([('waterway','dock')],5,'name'),
            ([],10,'icao'),
            ([],10,'iata'),  
            ([],10,'access'),
            ([],10,'name'),
            ([],10,'name:en'),
            ([],10,'sport'),
            ([],10,'intermittent'),
            ([],10,'seasonal'),
            ([],12,'bridge'),
            ([],12,'tunnel'),
            ([],13,'*')])
        
        
        
        extra_tags['highway'] = (1,[
            ([],6,'highway'),
            ([],6,'railway'),
            ([],6,'aerialway'),
            ([], 6, 'z_order'),
            ([], 10, 'bridge'),
            ([], 10, 'tunnel'),
            ([], 10, 'access'),
            ([('highway','motorway')], 4, 'ref'),
            ([('highway','trunk')], 4, 'ref'),
            ([('highway','primary')], 6, 'ref'),
            ([('highway','secondary')], 11, 'ref'),
            
            ([], 11, 'highspeed'),
            ([], 12, 'bus_routes'),
            ([], 12, 'cycle_routes'),
            ([], 12, 'service'),
            
            ([],13, '*'),
            
        ])
    
    extra_tags['polygon_exterior']=(1,extra_tags['polygon'][1])
    
    
    mzs = dict(((a,b,c),(d,[x for x in e.split(";") if x in extra_tags])) for a,b,c,d,e in oqt.geometry.minzoomvalues.default)
    
    if polypoint:
        
        extra_tags['polypoint']=(0,extra_tags['polygon'][1])
    
    for k,(va,vb) in mzs.items():
        if polypoint:
            if 'boundary' in vb or 'polygon' in vb:
                vb.append('polypoint')
                vb.append('polygon_exterior')
                mzs[k]=(va,vb)
        
        if not vb:
            vb.append('point' if a==0 else 'line' if a==1 else 'polygon')
            mzs[k]=(va,vb)
    
       
    
    return mzs, extra_tags
        
def prep_pa(poly, minzoom, cb,use_nt,polypoint,mergegeoms=True,alltags=False, simp_max=False,otherfeatures=None):
    feats,extra_tags = prep_spec(polypoint,alltags)
    
    if use_nt:
        return _oqttiles.make_processall_alt_callback_nt(gs, feats, extra_tags, minzoom,True,poly, True, mergegeoms, simp_max, otherfeatures,cb)
    return _oqttiles.make_processall_alt_callback(gs, feats, extra_tags, minzoom,True,poly, True, mergegeoms, simp_max, otherfeatures,cb)
    
def prep_mtd(filter_box,polypoint=False,alltags=False,simp_max=False):
    feats,extra_tags = prep_spec(polypoint,alltags)
    return _oqttiles.MakeTileData(feats, extra_tags, 14, True, True, filter_box, True, simp_max)




class testpp:
    def __init__(self, x,y,z,maxzoom):
        self.bounds=None
        if not x is None:
            q=oqt.elements.quadtree_from_tuple(x,y,z)
            self.bounds=oqt.elements.quadtree_bbox(q,0.05)
        
        self.maxzoom=maxzoom
    def __call__(self, q):
        if not self.maxzoom is None:
            if (q&31)>self.maxzoom:
                return False
        if self.bounds is None:
            return True
        return self.bounds.overlaps_quadtree(q)


MERGEGEOMS=True#False
ALLTAGS=False#True

class GeometryFromFile:
    def __init__(self, geomsfn):
        self.geomsfn=geomsfn
        
        self.hh=oqt.pbfformat.get_header_block(geomsfn).Index
    def __call__(self, x, y, z,maxzoom=None):
        return iter_geoms(self.geomsfn, self.hh, testpp(x,y,z,maxzoom),maxzoom)

class GeometryFromOqt:
    def __init__(self, prfx, timestamp=None):
        self.prfx = prfx
        self.timestamp=timestamp
    
    def __call__(self, x, y, z, maxzoom=None):
        bx = None if x is None else oqt.elements.quadtree_bbox(oqt.elements.quadtree_from_tuple(x,y,z),0.1)
        errs,geoms = oqt.process_geometry(self.prfx, bx, lastdate=self.timestamp, mergetiles=True, maxtilelevel=maxzoom)
        
        return [geoms]
        
class WriteToMbTiles:
    def __init__(self, fn, x,y,z, minzoom=None, maxzoom=None, timestamp=None, remove_existing=False):
        self.fn = fn
        self.minzoom=minzoom
        self.maxzoom=maxzoom
        
        self.total=0
        self.basetile = 'lowzoom' if x is None else tuple_to_str((x,y,z))
        self.timestamp=timestamp
        
        if remove_existing:
            tiles=mbt.MBTiles(self.fn)
            tiles.remove_basetile(self.basetile)
        
    def write_summary(self):
        if self.total>0:
            print("write summary %s %s %d %d %d" % (self.basetile, self.timestamp, self.minzoom, self.maxzoom, self.total))
            tiles=mbt.MBTiles(self.fn)
            tiles.start()
            tiles.write_summary(self.basetile, self.timestamp, self.minzoom, self.maxzoom, self.total)
            tiles.finish()
        
        
    def __call__(self, tls):
        if not tls:
            print('nothing to write')
            return
        print("write %d tiles to %s [%s]" % (len(tls),self.fn,self.basetile))
        tiles=mbt.MBTiles(self.fn)
        tiles.start()
        for (x,y,z),data in tls.iteritems():
            if (not self.minzoom is None) and z<self.minzoom:
                continue
            if (not self.maxzoom is None) and z>self.maxzoom:
                continue
            tile=tuple_to_str((x,y,z))
            tiles.add_tile(x,y,z,data,self.basetile,tile)
            self.total+=1
        tiles.finish()


class PrepTiles:
    def __init__(self, src, dbfn, timestamp=None, polypoint=True, mergegeoms=True, alltags=False):
        if not os.path.exists(dbfn):
            tiles=mbt.MBTiles(dbfn,create=True,minzoom=0,maxzoom=14,alt_schema=True)
        
        self.dbfn=dbfn
        
        if os.path.isdir(src):
            self.geoms = GeometryFromOqt(src,timestamp)
        else:
            self.geoms = GeometryFromFile(src)
        
        self.timestamp=timestamp
        self.polypoint=polypoint
        self.mergegeoms=mergegeoms
        self.alltags=alltags
        
    
    def write_lowzoom(self, maxzoom=9, otherfeatures=None, nothread=False):
        self(None,None,None,0,maxzoom, simp_max=True,otherfeatures=otherfeatures,nothread=nothread)
    
    
    def write_with_objs(self, objs, x, y, z, minzoom, maxzoom, simp_max, otherfeatures,nothread):
        st=time.time()
        tilepoly=box(*_oqttiles.tile_bound(0,0,0,0)) if x is None else box(*_oqttiles.tile_bound(x,y,z,-0.000001))
        
        tiles = WriteToMbTiles(self.dbfn, x, y, z, minzoom, maxzoom or 14, self.timestamp, True)
        pa=prep_pa(tilepoly, maxzoom or 14, tiles,nothread,self.polypoint,self.mergegeoms,self.alltags,simp_max,otherfeatures)
        
        nobjs=0
        for gg in objs:
            if not gg:
                continue
            nobjs+=sum(len(g) for g in gg)
            minq=min((g.Quadtree for g in gg),key=lambda q: q&31)
            print("[%9.1fs %6.1f%%] %-18s {%s => %s} %d tiles, %8d objs" % (time.time()-st,gg[-1].FileProgress,oqt.elements.quadtree_string(minq),"(%5d %5d %2d)" % oqt.elements.quadtree_tuple(gg[0].Quadtree), "(%5d %5d %2d)" % (oqt.elements.quadtree_tuple(gg[-1].Quadtree)), len(gg),nobjs))
            pa(gg)
        pa([])
        
        tiles.write_summary()
        print("[%7.1fs] %s: %d objs, %d tiles" % (time.time()-st,tiles.basetile,nobjs,tiles.total))
        return tiles.total
    
    
    def __call__(self, x, y, z, minzoom=10, maxzoom=None, simp_max=False, otherfeatures=None, nothread=False):
        
        return self.write_with_objs(self.geoms(x,y,z,maxzoom), x,y,z,minzoom,maxzoom,simp_max,otherfeatures,nothread)
        
        
        
        
        

def prep_tiles(geoms, x,y,z,callback, minzoom=8,maxzoom=None,use_nt=False,polypoint=False, mergegeoms=False, alltags=True):
    
    st=time.time()
    tilepoly=box(*_oqttiles.tile_bound(x,y,z,-0.000001))
    print("%d %d %d {%s}: %0.1fkm2 [%.80s]" % (x,y,z,tuple_to_str((x,y,z)),tilepoly.area/1000000,tilepoly.wkt))
    
    
    pa=prep_pa(tilepoly, maxzoom or 14, callback,use_nt,polypoint,mergegeoms,alltags)
    
    nobjs=0
    
    
    for gg in geoms:
        nobjs+=sum(len(g) for g in gg)
        minq=min((g.Quadtree for g in gg),key=lambda q: q&31)
        print("[%9.1fs %6.1f%%] %-18s {%s => %s} %d tiles, %8d objs" % (time.time()-st,gg[-1].FileProgress,oqt.elements.quadtree_string(minq),"(%5d %5d %2d)" % oqt.elements.quadtree_tuple(gg[0].Quadtree), "(%5d %5d %2d)" % (oqt.elements.quadtree_tuple(gg[-1].Quadtree)), len(gg),nobjs))
        pa(gg)
    pa([])
    return nobjs

def prep_tiles_lowzoom(geomsfn, hh, callback, maxzoom=7,use_nt=False,polypoint=False, mergegeoms=MERGEGEOMS, alltags=ALLTAGS):
    
    st=time.time()
    tilepoly=box(*_oqttiles.tile_bound(0,0,0,0))
    print("lowzoom {%s}: %0.1fkm2 [%.80s]" % (tuple_to_str((0,0,0)),tilepoly.area/1000000,tilepoly.wkt))
    
    
    pa=prep_pa(tilepoly, maxzoom or 14, callback,use_nt,polypoint,mergegeoms,alltags)
    
    nobjs=0
    for gg in iter_geoms(geomsfn, hh, None,maxzoom,25):
        nobjs+=sum(len(g) for g in gg)
        minq=min((g.Quadtree for g in gg),key=lambda q: q&31)
        print("[%9.1fs %6.1f%%] %-18s {%s => %s} %d tiles, %8d objs" % (time.time()-st,gg[-1].FileProgress,oqt.elements.quadtree_string(minq),"(%5d %5d %2d)" % oqt.elements.quadtree_tuple(gg[0].Quadtree), "(%5d %5d %2d)" % (oqt.elements.quadtree_tuple(gg[-1].Quadtree)), len(gg),nobjs))
        pa(gg)
    pa([])
    return nobjs
