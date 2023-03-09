def make_DF_vtk_cell(SKIP=1):
    import os
    import numpy as np
    import io
    import glob
    #import sys
    #
    # dt_out=sys.argv[1]
    # dt_out=int(sys.argv[1])
    # print('dt_out = ',dt_out)
    #dt_out = 600
    # SKIP = 6
    #
    ls='\n'
    #
    tpath = "output_DF_vtk"
    file_list = glob.glob("*")
    if tpath in file_list:
        tpath = os.path.join(os.getcwd(), tpath)
    else:
        os.mkdir(os.path.join(os.getcwd(), tpath))
        tpath = os.path.join(os.getcwd(), tpath)

    #srcpath=os.path.join('..','..','04_EXE_SIMU','02_DF_small_shift')
    srcpath = os.getcwd()
    varis=('h','cc','dz')
    fo_h='DF_{}_{:010d}.vtk'
    fz=os.path.join(srcpath,'elevation_inDF.asc')
    fb=os.path.join(srcpath,'targetArea_inDF.asc')
    files0=glob.glob(os.path.join(srcpath,'output_DF','res_*_??????????.asc'))
    files0.sort()
    #
    tmp_time1 = int(os.path.basename(files0[0])[-14:-4])
    tmp_time2 = int(os.path.basename(files0[1])[-14:-4])
    delta_time = tmp_time2 - tmp_time1
    dt_out = delta_time * SKIP
    #
    z=np.loadtxt(fz,skiprows=6,dtype=float)
    basin=np.loadtxt(fb,skiprows=6,dtype=int)
    for var in varis:
        files=[fn for fn in files0 if '_%s_'%var in os.path.basename(fn)]
        for fn in files:
            sec=int(os.path.basename(fn).split('_')[-1][-14:-4])
            if sec%dt_out==0:
                #fvtk=os.path.basename(fn)[:-4]+'.vtk'
                fvtk=fo_h.format(var,sec)
                fvtk=os.path.join(tpath ,fvtk)
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

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-s", "--skip", \
                    type=int, \
                    default=1, \
                    help = "read files in 'output_DR' by each SKIP interval"
                    )
parser.add_argument("--version", version="%(prog)s 1.0.0",
                    action="version",
                    default=False)
args = parser.parse_args()

if args.skip > 0:
    make_DF_vtk_cell(SKIP=args.skip)
else:
    import sys
    parser.print_help(sys.stderr)
    sys.exit(1)
