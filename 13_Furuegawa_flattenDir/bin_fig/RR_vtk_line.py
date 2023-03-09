def make_RR_vtk_line(SKIP=1):
    import os
    import numpy as np
    import io
    import glob
    # import sys
    #
    # dt_out=sys.argv[1]
    # dt_out=int(sys.argv[1])
    # print('dt_out = ',dt_out)
    #dt_out = 600
    # SKIP = 1
    #
    ls='\n'
    #
    tpath = "output_RR_vtk"
    file_list = glob.glob("*")
    if tpath in file_list:
        tpath = os.path.join(os.getcwd(), tpath)
    else:
        os.mkdir(os.path.join(os.getcwd(), tpath))
        tpath = os.path.join(os.getcwd(), tpath)

    # srcpath=os.path.join('..','..','04_EXE_SIMU','01_RR_DR')
    srcpath = os.getcwd()
    fo_h='RR_line_{0:010d}.vtk'
    files=glob.glob(os.path.join(srcpath,'output_RR','res_1d_??????????.txt'))
    files.sort()
    #
    tmp_time1 = int(os.path.basename(files[0])[-14:-4])
    tmp_time2 = int(os.path.basename(files[1])[-14:-4])
    delta_time = tmp_time2 - tmp_time1
    dt_out = delta_time * SKIP
    #
    fvtk=open('stream_line.vtk','w')
    fn=os.path.join(srcpath,'streamConfiguration_inRR.txt')
    flgs=np.loadtxt(fn,skiprows=2,delimiter=',',dtype=str)
    xa=flgs.T[12].astype(float)
    ya=flgs.T[13].astype(float)
    za=flgs.T[ 5].astype(float)
    ba=flgs.T[ 9].astype(float)
    acc=flgs.T[ 11].astype(int)
    acc=np.log10(acc)
    na=flgs.T[ 0].astype(int)
    npo=len(na)
    head=['# vtk DataFile Version 2.0',
        'Line example','ASCII','DATASET POLYDATA','POINTS %d float'%npo]
    io_head = io.StringIO()
    io_head.write(ls.join(head)+ls)
    for x,y,z in zip(xa,ya,za) :
        io_head.write('%f %f %f%s'%(x,y,z,ls))
    na2=flgs.T[1].astype(int)
    nli=np.size(np.where(na2!=0))
    io_head.write('LINES %d %d%s'%(nli,3*nli,ls))
    for i,i2 in zip(na,na2):
        if 0<i2:
            io_head.write('%d %d %d%s'%(2,i-1,i2-1,ls))
    head=['CELL_DATA %d'%nli,'SCALARS %s %s 1'%('z','float'),
        'LOOKUP_TABLE default']
    io_head.write(ls.join(head)+ls)
    itr=1
    for i,i2 in zip(na,na2):
        if 0<i2:
            z=za[i-1]
            if itr%10==0:
                io_head.write('%f%s'%(z,ls))
            else:
                io_head.write('%f '%(z))
        itr=itr+1
    for fn in files:
        sec=int(os.path.basename(fn).split('_')[-1][-14:-4])
        if sec%dt_out==0:
            data=np.loadtxt(fn,skiprows=0,dtype=str)
            za=data.T[3].astype(float)
            ha=data.T[4].astype(float)
            ua=data.T[5].astype(float)
            qa=data.T[6].astype(float)
            fvtk=fo_h.format(sec)
            fvtk=os.path.join(tpath,fvtk)
            print('output --->',fvtk)
            vtk=open(fvtk,'w')
            vtk.write(io_head.getvalue())
            for cvar,var in [['h',ha],['u',ua],['q',qa],['z',za]]:
                vtk.write(ls)
                head=['SCALARS %s %s 1'%(cvar,'float'),'LOOKUP_TABLE default']
                vtk.write(ls.join(head)+ls)
                itr=1
                for i,i2 in zip(na,na2):
                    if 0<i2:
                        v=var[i-1]
                        if itr%10==0:
                            vtk.write('%f%s'%(v,ls))
                        else:
                            vtk.write('%f '%(v))
                    itr=itr+1
            vtk.close()
    io_head.close()

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
    make_RR_vtk_line(SKIP=args.skip)
else:
    import sys
    parser.print_help(sys.stderr)
    sys.exit(1)
