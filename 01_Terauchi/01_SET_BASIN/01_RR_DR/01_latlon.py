import sys
import os
import glob
files=glob.glob(os.path.join('..','data_from_GSI','FG-GML-????-??-DEM10B.tif'))
res=10
print('input ---> ',sys.argv)
fdem='dem_%dm'%res
if len(sys.argv)==1:
    fns=[]
    for fn0 in files:
        fn=os.path.basename(fn0)
        print (fn)
        cmd='r.in.gdal -o in=%s out=%s --o'%(fn0,fn[:-4])
        os.system(cmd)
        fns.append(fn[:-4])
    cmd='g.region rast=%s'%(','.join(fns))
    os.system(cmd)
    cmd='r.patch in=%s out=%s --o'%(','.join(fns),fdem)
    os.system(cmd)
    cmd='r.watershed -a ele=%s drain=dir acc=acc --o'%fdem
    os.system(cmd)
    cmd='r.out.gdal in=acc out=acc_ll.tif --o'
    os.system(cmd)
    print('Step 1, output ---> acc_ll.tif')
elif len(sys.argv)==3:
    lon=sys.argv[1]
    lat=sys.argv[2]
    try:
        x=float(lon)
        y=float(lat)
        cmd='r.water.outlet input=dir output=basin coordinates=%f,%f --o'%(x,y)
        os.system(cmd)
        cmd='r.out.gdal in=basin out=basin_ll.tif --o'
        os.system(cmd)
        cmd='r.to.vect in=basin out=basin type=area --o'
        os.system(cmd)
        cmd='v.buffer in=basin out=basin_buf dis=%f --o'%(3./3600)# map unit
        os.system(cmd)
        print('Step 2, set outlet <--- (x,y)=({}, {})'.format(x,y))
        print('Step 2, output ---> basin_ll.tif')
    except:
        print('Format Error@(x,y) ---> exit')
        exit()
else:
    print('Error@input? ---> exit')




