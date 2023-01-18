import os
import numpy as np
import io
import glob
import sys
#
dt_out=sys.argv[1]
dt_out=int(sys.argv[1])
print('dt_out = ',dt_out)
#
ls='\n'
tpath=os.path.join('..','..','04_EXE_SIMU','02_DF_small_shift')
varis=('h','cc','dz')
fo_h='DF_{}_{:010d}.vtk'
fz=os.path.join(tpath,'input','elevation_inDF.asc')
fb=os.path.join(tpath,'input','targetArea_inDF.asc')
files0=glob.glob(os.path.join(tpath,'output_DF','res_*_??????????.asc'))
#
z=np.loadtxt(fz,skiprows=6,dtype=float)
basin=np.loadtxt(fb,skiprows=6,dtype=int)
for var in varis:
    files=[fn for fn in files0 if '_%s_'%var in os.path.basename(fn)]
    for fn in files:
        sec=int(os.path.basename(fn).split('_')[-1][4:-4])
        if sec%dt_out==0:
            #fvtk=os.path.basename(fn)[:-4]+'.vtk'
            fvtk=fo_h.format(var,sec)
            fvtk=os.path.join('out',fvtk)
            print('output --->',fvtk)
            vtk=open(fvtk,'w')
            head=open(fn).read().split(ls)[:6]
            ncols,nrows=int(head[0].split()[1]),int(head[1].split()[1])
            xco,yco=float(head[2].split()[1]),float(head[3].split()[1])
            cs=float(head[4].split()[1])
            #
            h=np.loadtxt(fn,skiprows=6,dtype=float)
            head=['# vtk DataFile Version 3.0','cell',
                  'ASCII','DATASET STRUCTURED_POINTS',
                  'DIMENSIONS %d %d %d'%(ncols,nrows,1),
                  'SPACING %f %f %f'%(cs,cs,0),
                  #'ORIGIN %f %f %f'%(xco+0.5*cs,yco+0.5*cs,0),
                  'ORIGIN %f %f %f'%(xco,yco,0),
                  'POINT_DATA %d'%(ncols*nrows),
                  'SCALARS %s %s 1'%('z','float'),
                  'LOOKUP_TABLE default']
            vtk.write(ls.join(head)+ls)
            io1 = io.StringIO()
            np.savetxt(io1,z[::-1])
            vtk.write(io1.getvalue())
            io1.close()
            #
            subhead=['SCALARS %s %s 1'%('basin','int'),'LOOKUP_TABLE default']
            vtk.write(ls.join(subhead)+ls)
            io1 = io.StringIO()
            np.savetxt(io1,basin[::-1],fmt='%d',delimiter=' ')
            vtk.write(io1.getvalue())
            #
            subhead=['SCALARS %s %s 1'%(var,'float'),'LOOKUP_TABLE default']
            vtk.write(ls.join(subhead)+ls)
            h=np.loadtxt(fn,skiprows=6,dtype=float)
            io2 = io.StringIO()
            np.savetxt(io2,h[::-1],fmt='%f',delimiter=' ')
            vtk.write(io2.getvalue())
            #
            vtk.close()
