def make_DF_minmaxasc():
    import os
    import glob
    import sys
    import numpy as np
    file_list = glob.glob("*")
    if "output_DF_asc" in file_list:
        tdir = os.path.join(os.getcwd(), "output_DF_asc")
    else:
        os.mkdir(os.path.join(os.getcwd(), "output_DF_asc"))
        tdir = os.path.join(os.getcwd(), "output_DF_asc")
    print("DF out @ [dz, cc]")
    ls='\n'
    sdir = os.path.join(os.getcwd(), "output_DF")
    files_dz=glob.glob(os.path.join(sdir,'res_dz_??????????.asc'))
    files_dz.sort()
    files_cc=glob.glob(os.path.join(sdir,'res_cc_??????????.asc'))
    files_cc.sort()
    #
    fn=files_dz[0]
    head=open(fn).read().split(ls)[:6]
    ncols=int(head[0].split()[1])
    nrows=int(head[1].split()[1])
    xllcorner=float(head[2].split()[1])
    yllcorner=float(head[3].split()[1])
    dx=float(head[4].split()[1])
    dy=dx
    head[5] = head[5].replace("NODATA_value -9999.000", "NODATA_value 0.0")
    head=ls.join(head)
    #
    dzmin = np.zeros((nrows,ncols))
    dzmax = np.zeros((nrows,ncols))
    ccmax = np.zeros((nrows,ncols))
    #
    for i in range(len(files_dz)):
        dz = np.loadtxt(files_dz[i],skiprows=6)
        dzmax = np.maximum(dzmax, dz)
        dzmin = np.minimum(dzmin, dz)
        cc = np.loadtxt(files_cc[i],skiprows=6)
        ccmax = np.maximum(ccmax, cc)
    #
    np.savetxt(os.path.join(tdir,'dzmax_inDF.asc'),dzmax,fmt='%7.2f',header=head,comments='')
    np.savetxt(os.path.join(tdir,'dzmin_inDF.asc'),dzmin,fmt='%7.2f',header=head,comments='')
    np.savetxt(os.path.join(tdir,'ccmax_inDF.asc'),ccmax,fmt='%7.2f',header=head,comments='')
    #
    # rearrange .asc files in output_DF 
    # revise NODATA_value and copy
    names = ["rain_map.asc", "rain_sum.asc", "res_hmax.asc", "res_hsmax.asc"]
    for name in names:
        data = np.loadtxt(os.path.join(sdir, name), skiprows=6)
        if name == "rain_map.asc":
            name = "mask.asc"
        np.savetxt(os.path.join(tdir,name),data,fmt='%7.2f',header=head,comments='')


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--version", version="%(prog)s 1.0.0",
                    action="version",
                    default=False)
args = parser.parse_args()

try:
    make_DF_minmaxasc()
except Exception as e:
    import sys, traceback
    print(traceback.format_exc())
    sys.exit(-1)
