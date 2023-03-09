def make_DR_maxasc(SKIP=1):
    import os
    import numpy as np
    import glob
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    # DT should be set as the same value of output_f_time in DR_input.txt
    #DT = 3600
    #SKIP = 36

    file_list = glob.glob("*")
    # if "output_DR_png" in file_list:
    #     pass
    # else:
    #     os.mkdir(os.path.join(os.getcwd(), "output_DR_png"))
    if "output_DR_asc" in file_list:
        pass
    else:
        os.mkdir(os.path.join(os.getcwd(), "output_DR_asc"))
    #
    ls='\n'
    tdir=os.getcwd()
    files=glob.glob(os.path.join(tdir,'output_DR','res_??????????.txt'))
    files.sort()
    #
    tmp_time1 = int(os.path.basename(files[0])[-14:-4])
    tmp_time2 = int(os.path.basename(files[1])[-14:-4])
    delta_time = tmp_time2 - tmp_time1
    DT = delta_time * SKIP
    #
    fn=os.path.join(tdir,'targetArea_inRR.asc')
    head=open(fn).read().split(ls)[:6]
    xllcorner=float(head[2].split()[1])
    yllcorner=float(head[3].split()[1])
    dx=float(head[4].split()[1])
    dy=dx
    head=ls.join(head)
    #
    ij2d=np.loadtxt(fn,skiprows=6)
    ij2d=np.full(np.shape(ij2d),0,dtype=int)
    jend,iend=np.shape(ij2d)
    x=np.linspace(xllcorner,xllcorner+iend*dx,iend)
    y=np.linspace(yllcorner,yllcorner+jend*dy,jend)
    Xg,Yg=np.meshgrid(x,y)
    fn=os.path.join(tdir,'streamConfiguration_inRR.txt')
    idxs=np.loadtxt(fn,skiprows=2,delimiter=',',dtype=int,usecols=(14,15))
    i_id=np.transpose(idxs)[0]
    j_id=np.transpose(idxs)[1]
    #
    a=np.full(np.shape(ij2d),0,dtype=float)
    qmax=np.full(np.shape(ij2d),0,dtype=float)
    umax=np.full(np.shape(ij2d),0,dtype=float)
    hmax=np.full(np.shape(ij2d),0,dtype=float)
    dze=np.full(np.shape(ij2d),0,dtype=float)
    ifld=np.full(np.shape(ij2d),0,dtype=int)
    hsc=np.full(np.shape(ij2d),0,dtype=float)
    for fn in files[0:]:
        print (fn,)
        ucols=tuple(range(17))
        data=np.loadtxt(fn,skiprows=1,usecols=ucols)
        za=np.transpose(data)[1]
        ha=np.transpose(data)[2]
        ua=np.transpose(data)[3]
        zinia=np.transpose(data)[9]
        dza=za-zinia
        iflda=np.transpose(data)[15]
        hsca=np.transpose(data)[14]
        qa=np.transpose(data)[4]
        dero=np.sum(np.where(dza<0,dza,0))
        ddep=np.sum(np.where(dza>0,dza,0))
        #print ('%10.1f%10.1f%10.1f'%(dero,ddep,ddep-dero))
        #for i,j,h,u,dz,hs,qs,ifld0,hsc0,q in zip(i_id,j_id,ha,ua,dza,hsa,qsa,iflda,hsca,qa):
        for i,j,h,u,dz,ifld0,hsc0,q in zip(i_id,j_id,ha,ua,dza,iflda,hsca,qa):
            if i==0:
                print (i,j)
                continue
            umax[jend-j,i-1]=max(u,umax[jend-j,i-1])
            hmax[jend-j,i-1]=max(h,hmax[jend-j,i-1])
            dze[jend-j,i-1]=dz
            ifld[jend-j,i-1]=max(ifld0,ifld[jend-j,i-1])
            hsc[jend-j,i-1]=max(hsc0,hsc[jend-j,i-1])
            qmax[jend-j,i-1]=max(q,qmax[jend-j,i-1])
        # isec=int(os.path.basename(fn)[-14:-4])
        # if isec%DT==0:
        #     levels = [-1,0,1]
        #     fign=os.path.join('output_DR_png','dz_%s.png'%(os.path.basename(fn)[-14:-4]))
        #     plt.pcolor(Xg,Yg[::-1],dze,cmap='coolwarm', norm=Normalize(vmin=-1, vmax=1), shading='auto')
        #     plt.xlabel('East-West (m)')
        #     plt.ylabel('Sauth-North (m)')
        #     plt.colorbar(label='dz (m)')
        #     # plt.axes().set_aspect(dx/dy)
        #     plt.savefig(fign)
        #     plt.clf()
        #     plt.close()
        #     fign=os.path.join('output_DR_asc','dz_%s.asc'%(os.path.basename(fn)[-14:-4]))
        #     np.savetxt(fign,dze,fmt='%7.3f',header=head,comments='')
    np.savetxt(os.path.join("output_DR_asc",'umax.asc'),umax,fmt='%f',header=head,comments='')
    np.savetxt(os.path.join("output_DR_asc",'hmax.asc'),hmax,fmt='%7.2f',header=head,comments='')
    np.savetxt(os.path.join("output_DR_asc",'dz.asc'),dze,fmt='%6.2f',header=head,comments='')
    np.savetxt(os.path.join("output_DR_asc",'qmax.asc'), qmax, fmt='%f',header=head,comments='')
    #raw_input()

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
    make_DR_maxasc(SKIP=args.skip)
else:
    import sys
    parser.print_help(sys.stderr)
    sys.exit(1)

