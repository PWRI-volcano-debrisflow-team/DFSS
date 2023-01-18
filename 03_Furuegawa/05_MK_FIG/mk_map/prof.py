import os
import numpy as np
import glob
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

ls='\n'

tdir=os.path.join('..','..','04_EXE_SIMU','01_RR_DR')

fn=os.path.join(tdir,'input','targetArea_inRR.asc')
head=open(fn).read().split(ls)[:6]

xllcorner=float(head[2].split()[1])
yllcorner=float(head[3].split()[1])
dx=float(head[4].split()[1])
dy=dx

head=ls.join(head)
ij2d=np.loadtxt(fn,skiprows=6)
ij2d=np.full(np.shape(ij2d),0,dtype=np.int)
jend,iend=np.shape(ij2d)

x=np.linspace(xllcorner,xllcorner+iend*dx,iend)
y=np.linspace(yllcorner,yllcorner+jend*dy,jend)
Xg,Yg=np.meshgrid(x,y)

fn=os.path.join(tdir,'input','streamConfiguration_inRR.txt')
idxs=np.loadtxt(fn,skiprows=2,delimiter=',',dtype=np.int,usecols=(14,15))
i_id=np.transpose(idxs)[0]
j_id=np.transpose(idxs)[1]

idxs=np.loadtxt(fn,skiprows=2,delimiter=',',dtype=np.int,usecols=(0,1))
id_1d=np.transpose(idxs)[0]
ito_1d=np.transpose(idxs)[1]
idxs=np.loadtxt(fn,skiprows=2,delimiter=',',usecols=(7))
dx_1d=np.transpose(idxs)

idxs=np.loadtxt(fn,skiprows=2,delimiter=',',usecols=(7))
dx_1d=np.transpose(idxs)

files=glob.glob(os.path.join(tdir,'output_DR','res_??????????.txt'))
files.sort()

#files=[afile for afile in files if int(os.path.basename(afile)[4:-4])%600==0]

files=[afile for afile in files if int(os.path.basename(afile)[4:-4])%3600==0]

i_sta=8
i_end=1574

i0=i_sta
idxs=[i0]

while i0 != i_end:
    i1=ito_1d[np.where(id_1d==i0)][0]
    i0=i1
    idxs.append(i1)

dis0=[dx_1d[i-1] for i in idxs]
dis0=dis0[:]
dis=[sum(dis0[1:i+1]) for i in range(len(dis0))]

iofs=0
for if0,fn in enumerate(files[iofs:]):
    print (if0+iofs,fn)
    sec=int(os.path.basename(fn)[4:-4])
    data=np.loadtxt(fn,skiprows=1)
    za=np.transpose(data)[1]
    ha=np.transpose(data)[2]
    ua=np.transpose(data)[3]
    zinia=np.transpose(data)[9]
    dza=za-zinia
    hsa=np.transpose(data)[11]
    qsa=np.transpose(data)[12]
    iflda=np.transpose(data)[15]
    hsca=np.transpose(data)[14]
    qa=np.transpose(data)[16]
    
    z1d=[za[i-1] for i in idxs]
    z01d=[zinia[i-1] for i in idxs]
    hs1d=[(zinia+hsa)[i-1] for i in idxs]
    h1d=[(zinia+hsa+ha)[i-1] for i in idxs]    
        
    plt.plot(dis,z1d,lw=0.5,label='z {}s'.format(sec))
    #plt.plot(dis,hs1d,lw=1.,label='hs')
    #plt.plot(dis,h1d,lw=1.,label='h')    
    #plt.title('No.%d: %s'%(if0,os.path.basename(fn)))
    #if if0==0: plt.plot(z01d)
    #plt.xlim(1000,1500)
    #plt.ylim(100,170)
   
plt.plot(dis,z01d,lw=1,color='k',label='zini')
plt.xlabel('Distance (m)')
plt.ylabel('Elevation (m)')
plt.legend()
plt.show()




#dz=np.loadtxt('dz.asc',skiprows=6)
#dz_1d=[]
#for i,j in zip(i_id,j_id):
#    dz_1d.append(dz[jend-j,i-1])
#aaa=np.array(dz_1d)
#print np.nansum(np.where(aaa<0.,aaa,np.nan))




