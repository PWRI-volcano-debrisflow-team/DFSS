import sys
import os
import glob
def trs(lon,lat,src_EPSG,dst_EPSG):
    from osgeo import ogr, osr, gdal
    src_srs, dst_srs = osr.SpatialReference(), osr.SpatialReference()
    src_srs.ImportFromEPSG(src_EPSG)
    dst_srs.ImportFromEPSG(dst_EPSG)
    trans = osr.CoordinateTransformation(src_srs, dst_srs)
    return trans.TransformPoint(lat,lon)[::-1][1:3] # lat,lon -> lon,lat
res_ll=10
res=10
loc='JGD2011'
print('input ---> ',sys.argv)
fdem='dem_%dm'%res_ll
if len(sys.argv)==1:
    cmd='v.proj in=%s out=area loc=%s --o'%('basin_buf',loc)
    os.system(cmd)
    cmd='g.region -a vect=area res=%d'%(res)
    os.system(cmd)
    cmd='r.proj in=%s out=dem_%dm  loc=%s method=bilinear --o'%(fdem,res,loc)
    os.system(cmd)
    cmd='r.watershed -a ele=dem_%dm acc=acc drain=dir --o'%res
    os.system(cmd)
    cmd='r.fill.dir input=dem_%dm dir=dir output=dem_filled_%dm  --o'%(res,res)
    os.system(cmd)
    cmd='r.watershed -a ele=dem_filled_%dm acc=acc drain=dir --o'%res
    os.system(cmd)
    cmd='r.out.gdal in=acc out=acc_xy.tif --o'
    os.system(cmd)
    print('Step 1, output ---> acc_xy.tif')
elif len(sys.argv)==3:
    #x=sys.argv[1]
    #y=sys.argv[2]
    ### --->
    lon=sys.argv[1]
    lat=sys.argv[2]
    ### <---
    try:
        #x=float(x)
        #y=float(y)
        ### --->
        lon=float(lon)
        lat=float(lat)
        (x,y)=trs(lon,lat,6668,6670)
        print('({},{}) ---> ({},{})'.format(lon,lat,x,y))
        ### <---
        cmd='g.region -a vect=area res=%d'%(res)
        os.system(cmd)
        dem='dem_filled_%dm'%res
        cmd='r.water.outlet input=dir output=basin coordinates=%f,%f --o'%(x,y)
        os.system(cmd)
        cmd='r.buffer in=basin out=basin_buf dis=%d --o'%res
        os.system(cmd)
        cmd='g.region zoom=basin_buf'
        os.system(cmd)
        for f in [['dir','dir_%dm'%res,0],[dem,'dem_%dm'%res,-999]]:
            fi=f[0]
            fo=f[1]
            nv=f[2]
            os.system('r.out.gdal in=%s out=%s.asc format=AAIGrid nodata=%d --o'%(fi,fo,nv))
        print('Step 2, set outlet <--- (x,y)=({}, {})'.format(x,y))
        print('Step 2, output ---> dem_%dm.asc, dir_%dm.asc'%(res,res))
    except:
        print('Format Error@(x,y) ---> exit')
        exit()
else:
    print('Error@input? ---> exit')




