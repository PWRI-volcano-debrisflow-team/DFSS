import os
res_ll=10
res=10
area='area.shp'
loc='JGD2011'
fdem='dem_%dm'%res_ll
cmd='v.in.ogr in=%s out=area_DF --o'%(area)
os.system(cmd)
cmd='v.buffer in=area_DF out=area_buf dis=%d --o'%(res)
os.system(cmd)
cmd='g.region -a vect=area_buf res=%d'%(res)
os.system(cmd)
cmd='v.to.rast in=area_DF out=area_DF use=val val=1 --o'
os.system(cmd)
fdem2='dem_%dm'%res
cmd='r.proj in=%s out=%s loc=%s method=bilinear --o'%(fdem,fdem2,loc)
os.system(cmd)
cmd='r.watershed -a ele=%s acc=acc drain=dir --o'%(fdem2)
os.system(cmd)
for f in [['dir','dir_%dm'%res,0],[fdem2,'dem_%dm'%res,-999],['area_DF','basin_%dm'%res,0]]:
    fi=f[0]
    fo=f[1]
    nv=f[2]
    os.system('r.out.gdal in=%s out=%s.asc format=AAIGrid nodata=%d --o'%(fi,fo,nv))
exit()

