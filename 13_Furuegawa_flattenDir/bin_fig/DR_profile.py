def make_timeSeriesProfiles(start_id=2100, end_id=307, skip=6):
    import os
    import numpy as np
    import glob
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    # 全てを書き出すと時間解像度が細かすぎるため適当に調整する
    # DT = 3600
    # 上流側のriver_node id
    # i_sta=1574
    # 下流側のriver_node id
    # i_end=2406

    i_sta = start_id
    i_end = end_id
    SKIP = skip

    ls='\n'
    sdir=os.getcwd()
    ##################################
    tdir = "output_DR_fig"
    file_list = glob.glob("*")
    if tdir in file_list:
        tdir = os.path.join(os.getcwd(), tdir)
    else:
        tdir = os.path.join(os.getcwd(), tdir)
        os.mkdir(tdir)
    ##################################
    fn=os.path.join(sdir,'targetArea_inRR.asc')
    head=open(fn).read().split(ls)[:6]

    xllcorner=float(head[2].split()[1])
    yllcorner=float(head[3].split()[1])
    cellsize=float(head[4].split()[1])
    dx = cellsize
    dy=dx

    head=ls.join(head)
    ij2d=np.loadtxt(fn,skiprows=6)
    ij2d=np.full(np.shape(ij2d),0,dtype=np.int64) # zeros と同じ
    jend,iend=np.shape(ij2d)

    x=np.linspace(xllcorner,xllcorner+iend*dx,iend)
    y=np.linspace(yllcorner,yllcorner+jend*dy,jend)
    Xg,Yg=np.meshgrid(x,y)

    fn=os.path.join(sdir,'streamConfiguration_inRR.txt')
    # i_2d,j_2d
    idxs=np.loadtxt(fn,skiprows=2,delimiter=',',dtype=np.int64,usecols=(14,15))
    i_id=np.transpose(idxs)[0]
    j_id=np.transpose(idxs)[1]

    # id,id_to
    idxs=np.loadtxt(fn,skiprows=2,delimiter=',',dtype=np.int64,usecols=(0,1))
    id_1d=np.transpose(idxs)[0]
    ito_1d=np.transpose(idxs)[1]
    # zs_1d
    #idxs=np.loadtxt(fn,skiprows=2,delimiter=',',usecols=(7))
    idxs=np.loadtxt(fn,skiprows=2,delimiter=',',usecols=(7))
    dx_1d=np.transpose(idxs)

    idxs=np.loadtxt(fn,skiprows=2,delimiter=',',usecols=(7))
    dx_1d=np.transpose(idxs)

    files=glob.glob(os.path.join(sdir,'output_DR','res_??????????.txt'))
    files.sort()
    #
    tmp_time1 = int(os.path.basename(files[0])[-14:-4])
    tmp_time2 = int(os.path.basename(files[1])[-14:-4])
    delta_time = tmp_time2 - tmp_time1
    # print(tmp_time1, tmp_time2, delta_time)
    #files=[afile for afile in files if int(os.path.basename(afile)[-14:-4])%600==0]
    files=[afile for afile in files if int(os.path.basename(afile)[-14:-4])%(delta_time*SKIP)==0]

    i0=i_sta
    idxs=[i0]

    while i0 != i_end:
        i1=ito_1d[np.where(id_1d==i0)][0]
        i0=i1
        idxs.append(i1)

    dis0=[dx_1d[i-1] for i in idxs]
    dis0=dis0[:]
    dis=[sum(dis0[1:i+1]) for i in range(len(dis0))]

    time_sec = []
    z01d_array = []
    iofs=0
    for if0,fn in enumerate(files[iofs:]):
        print (if0+iofs,fn)
        sec=int(os.path.basename(fn)[-14:-4])
        time_sec.append(sec)
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
        z01d_array.append(z01d)

    stack_tmp = np.append(dis,z01d_array)
    nrows = np.int64(len(stack_tmp)/len(dis))
    stack = stack_tmp.reshape((nrows, len(dis)))
    fname = "DR_profile_fm_{}_to_{}.csv".format(str(i_sta), str(i_end))
    dataframe = np.transpose(stack)
    header_label = "x, z_sec=" + ", z_sec=".join( np.array(time_sec).astype(str) )
    format_label = "%.5f, " + ", ".join( ["%.5f" for i in range(len(time_sec))] )
    # print(header_label)
    # print(format_label)
    np.savetxt(os.path.join(tdir,fname), dataframe, delimiter=",", \
         header=header_label, fmt=format_label)
    #
    plt.plot(dis,z01d,lw=1,color='k',label='zini')
    plt.xlabel('Distance (m)')
    plt.ylabel('Elevation (m)')
    plt.legend()
    # plt.show()
    plt.savefig(os.path.join(tdir, 'DR_profile.png'), dpi=300)
    print("")
    print("make png file in {}".format(os.path.join(tdir, 'DR_profile.png')))

    #dz=np.loadtxt('dz.asc',skiprows=6)
    #dz_1d=[]
    #for i,j in zip(i_id,j_id):
    #    dz_1d.append(dz[jend-j,i-1])
    #aaa=np.array(dz_1d)
    #print np.nansum(np.where(aaa<0.,aaa,np.nan))


def make_csvfile(dis, value_array, label="z"):
    stack_tmp = np.append(dis,z01d_array)
    nrows = np.int64(len(stack_tmp)/len(dis))
    stack = stack_tmp.reshape((nrows, len(dis)))
    fname = "DR_profile_fm_{}_to_{}.csv".format(str(i_sta), str(i_end))
    dataframe = np.transpose(stack)
    header_label = "x, z_sec=" + ", z_sec=".join( np.array(time_sec).astype(str) )
    format_label = "%.5f, " + ", ".join( ["%.5f" for i in range(len(time_sec))] )
    print(header_label)
    print(format_label)
    np.savetxt(os.path.join(tdir,fname), dataframe, delimiter=",", \
         header=header_label, fmt=format_label)


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("start_node_id", \
                    type=int, \
                    default=0, \
                    help="river node id of starting point where in upper region; ex) 1574", \
                    )
parser.add_argument("end_node_id", \
                    type=int, \
                    default=0, \
                    help="river node id of ending point where in upper region; ex) 2406", \
                    )
parser.add_argument("-s", "--skip", \
                    type=int, \
                    default=1, \
                    help = "read files in 'output_DR' by each SKIP interval"
                    )
parser.add_argument("--version", version="%(prog)s 1.0.0",
                    action="version",
                    default=False)
args = parser.parse_args()

if args.start_node_id != 0 and args.end_node_id != 0:
    make_timeSeriesProfiles(start_id=args.start_node_id, end_id=args.end_node_id, skip=args.skip)
else:
    import sys
    parse.print_help(sys.stderr)
    sys.exit(1)
