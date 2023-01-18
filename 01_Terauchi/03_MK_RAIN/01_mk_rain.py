import zipfile
import io
import numpy as np
import sys
import datetime
import os
import glob
import matplotlib.pyplot as plt
#
def aoi(rasterfn):
    #https://pcjericks.github.io/py-gdalogr-cookbook/raster_layers.html#get-raster-band
    #import gdal
    from osgeo import gdal
    raster = gdal.Open(rasterfn)
    geotransform = raster.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]
    band = raster.GetRasterBand(1)
    array = band.ReadAsArray()
    cols = array.shape[1]
    rows = array.shape[0]
    xmin,xmax=originX,originX+pixelWidth*cols
    ymin,ymax=originY+pixelHeight*rows,originY
    return(xmin,xmax,ymin,ymax)
#
def eva(doy):
    import math
    e0,m,ofs,lam=2.6,1.7,205,365
    e=math.cos((doy-ofs)/lam*math.pi*2.)*m+e0
    e=e/24.
#
    e=0.
#
    return e
#
def trs(lon,lat,src_EPSG,dst_EPSG):
    from osgeo import ogr, osr, gdal
    src_srs, dst_srs = osr.SpatialReference(), osr.SpatialReference()
    src_srs.ImportFromEPSG(src_EPSG)
    dst_srs.ImportFromEPSG(dst_EPSG)
    trans = osr.CoordinateTransformation(src_srs, dst_srs)
    return trans.TransformPoint(lat,lon)[::-1][1:3] # lat,lon -> lon,lat
#
ROI_xmin,ROI_xmax,ROI_ymin,ROI_ymax=aoi('targetArea_inRR.asc')
print('ROI(xy)',ROI_xmin,ROI_xmax,ROI_ymin,ROI_ymax)
pnts = [(ROI_xmin, ROI_ymin), (ROI_xmax, ROI_ymax)]
pnts=[trs(*pnt,6670,6668) for pnt in pnts]
[(ROI_xmin, ROI_ymin), (ROI_xmax, ROI_ymax)]=pnts
print('ROI(ll)',ROI_xmin,ROI_xmax,ROI_ymin,ROI_ymax)
files=glob.glob(os.path.join('XRAIN_from_DIAS','201707040000-201707080000-10-FUK-130.7031-33.5000-130.8500-33.3688.zip'))
files.sort()
for fn in files:
    fn_ps='_'.join(os.path.basename(fn).split('-')[0:2])
    myzip=zipfile.ZipFile(fn, mode='r')
    fs=','
    ls='\n'
    dt_min=float(os.path.basename(fn).split('-')[2])
    print('dt=',dt_min)
    xmin,ymin=[float(fn[:-4].split('-')[i]) for i in [4,7]]
    xmax,ymax=[float(fn[:-4].split('-')[i]) for i in [6,5]]
    fname=myzip.namelist()[0]
    myfile=myzip.open(fname)
    zipcont = io.StringIO(myfile.read().decode())
    output=np.loadtxt(zipcont,dtype=np.float,delimiter=fs)
    jend,iend=np.shape(output)
    dx,dy=(xmax-xmin)/float(iend-1),(ymax-ymin)/float(jend-1)
    xmin=xmin+0.5*dx
    ymin=ymin+0.5*dy
    print('AOI(ll)',xmin,xmin+float(iend)*dx,ymin,ymin+dy*float(jend))
    #print(dx,dy)
    x=np.linspace(xmin,xmin+(iend-1)*dx,iend)
    y=np.linspace(ymin,ymin+(jend-1)*dy,jend)
    X,Y=np.meshgrid(x,y)
    ij=np.where((ROI_xmin-dx < X) & (X < ROI_xmax+dx) & \
                (ROI_ymin-dy < Y) & (Y < ROI_ymax+dy))
    data=(np.vstack((X[ij],Y[ij]))).T
    xy=np.array([trs(*pnt,6668,6670) for pnt in data])
    fo=open('grid.csv','w')
    fo.write('lon,lat,x,y%s'%ls)
    for (lat,lon),(x0,y0) in zip(data,xy):
        #print(lat,lon,x0,y0)
        fo.write('%f,%f,%f,%f%s'%(lat,lon,x0,y0,ls))
    fo.close()
    fo=open('rain_time_%s.txt'%fn_ps,'w')
    fo.write('%d%s'%(len(xy.T[0]),ls))
    fo.write('x,'+fs.join(['%f'%f for f in xy.T[0]])+ls)
    fo.write('y,'+fs.join(['%f'%f for f in xy.T[1]])+ls)
    with zipfile.ZipFile(fn) as myzip:
        for i0,fn in enumerate(myzip.namelist()[:]):
            year,month,day,hour,minu=[int(fn[i1:i2]) for i1,i2 in [[0,4],[4,6],[6,8],[9,11],[11,13]]]
            date=datetime.datetime(year,month,day,hour,minu)
            if i0==0: date0=date
            #print(fn,date)
            with myzip.open(fn) as myfile:
                zipcont = io.StringIO(myfile.read().decode())
                output = np.loadtxt(zipcont,dtype=np.float,delimiter=fs)
                output=output[::-1]
                if date == datetime.datetime(2017,7,5,15,10):
                    plt.pcolor(X,Y,np.where(output<0.,np.nan,output),cmap=plt.get_cmap('jet'))
                    plt.colorbar(label='Rain (mm)')
                    plt.contour(X,Y,np.where(output<0.,np.nan,output),levels=[50,100],colors=['white'])
                    plt.axes().set_aspect(dx/dy)
                    #plt.show()
                    plt.savefig('{0:%Y%m%d%H%M.png}'.format(date),dpi=300)
                R2=output[ij]
                if min(R2)<0.:
                    print('{} ---> '.format(min(R2)),end='')
                    R2=np.where(R2<0.,0.,R2)
                    print('{}   '.format(min(R2)),'{0:%Y/%m/%d %H:%M}'.format(date))
                R2=R2*dt_min/60. # mm/h -> mm/dt
                doy=(date-datetime.datetime(date.year,1,1)).days+hour/24.
                e=eva(doy)
                nl=fs.join(['%f'%f for f in R2])
                nl="{0:'%Y/%m/%d %H:%M',}".format(date)+nl+',%f'%e
                #print (nl)
                fo.write(nl+ls)
            if (datetime.timedelta(minutes=dt_min) < date-date0):
                print('Leap!',date0,'--->',date)
            date0=date
    fo.close()
