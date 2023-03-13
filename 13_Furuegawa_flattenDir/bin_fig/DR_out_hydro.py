# -*- coding: utf-8 -*-
# def make_DR_fig(nodes, start_date, end_date, skip):
def make_DR_fig(loc=[], SKIP=1):
    import os
    import glob
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import datetime
    from matplotlib.dates import DateFormatter
    from matplotlib.ticker import ScalarFormatter
    ##################################
    r_gold=(1.+5.**0.5)*0.5
    r_silver=(1.+2.**0.5)
    r_yamato=2.**0.5
    w_inchi=3.1496
    h_inchi=w_inchi
    r_top=0
    r_bottom=0
    r_left=0
    r_right=0
    fsize=7
    ##################################
    #plt.figure(figsize=(w_inchi,h_inchi))
    plt.rcParams["font.size"] = fsize
    # from matplotlib.font_manager import FontProperties
    # fp = FontProperties(fname=r'C:\WINDOWS\Fonts\msgothic.ttc')
    ##################################
    tdir = "output_DR_fig"
    file_list = glob.glob("*")
    if tdir in file_list:
        tdir = os.path.join(os.getcwd(), tdir)
    else:
        tdir = os.path.join(os.getcwd(), tdir)
        os.mkdir(tdir)
    ##################################
    # fig, ax = plt.subplots(3,1, figsize=(w_inchi,h_inchi))
    # for i in range(3):
    #     ax[i].xaxis.set_major_formatter(DateFormatter('%b.%d\n%H:%M'))
    # ax2=[ax[i].twinx() for i in range(3)]
    # [ax2ax.xaxis.set_major_formatter(DateFormatter('%b.%d\n%H:%M')) for ax2ax in ax2]
    ls='\n'
    sdirs=glob.glob(os.getcwd())
    # loc=sys.argv[1:]
    loc=[int(s) for s in loc]
    print('DR out @',loc)
    lgs=[[1] for i in range(len(sdirs))]
    for lg,sdir in zip(lgs,sdirs):
        files=glob.glob(os.path.join(sdir,'output_DR','res_*.txt'))
        files.sort()
        tmp_time1 = int(os.path.basename(files[0])[-14:-4])
        tmp_time2 = int(os.path.basename(files[1])[-14:-4])
        delta_time = tmp_time2 - tmp_time1
        files=[afile for afile in files if int(os.path.basename(afile)[-14:-4])%(delta_time*SKIP)==0]
        iend=len(files)
        label=sdir.split('_')[-1]
        fn=os.path.join(sdir,'output_RR_misc','basin_rain.txt')
        cont=open(fn).read().split(ls)[1].split(',')[0].strip("'")
        date,time=cont.split()[0].split('/'),cont.split()[1].split(':')
        year,month,day=int(date[0]),int(date[1]),int(date[2])
        hour,minute=int(time[0]),int(time[1])
        time0=datetime.datetime(year,month,day,hour,minute)
        time0=time0-datetime.timedelta(seconds=0)
        hydro,hhh,uuu=[],[],[]
        qdf,pcc,pcf=[],[],[]
        qa,qca,qfa=[],[],[]
        time=[]
        time_sec=[]
        for i0,afile in enumerate(files[:]):
            #print (afile)
            res=np.loadtxt(afile,usecols=(4,5,6,16,17,18),skiprows=1)
            sec=float(os.path.basename(afile).split('_')[-1][:-4])
            #
            q=np.transpose(res)[0]
            cc=np.transpose(res)[1]
            cf=np.transpose(res)[2]
            qdf.append([q[l-1] for l in loc])
            pcc.append([cc[l-1] for l in loc])
            pcf.append([cf[l-1] for l in loc])
            #
            qaa=np.transpose(res)[3]
            qcaa=np.transpose(res)[4]
            qfaa=np.transpose(res)[5]
            qa.append([qaa[l-1] for l in loc])
            qca.append([qcaa[l-1] for l in loc])
            qfa.append([qfaa[l-1] for l in loc])
            #
            time.append(time0+datetime.timedelta(seconds=sec))
            time_sec.append(sec)
        for i in range(len(loc)):
            dt=(time[1]-time[0]).seconds
            qa,qca,qfa=np.array(qa),np.array(qca),np.array(qfa)
            q=np.append([0],(qa.T[i][1:]-qa.T[i][:-1])/dt)
            qc=np.append([0],(qca.T[i][1:]-qca.T[i][:-1])/dt)
            qf=np.append([0],(qfa.T[i][1:]-qfa.T[i][:-1])/dt)
            #
            cc = np.where(q > 0.0, qc/q, 0.0)
            cf = np.where(q > 0.0, qf/q, 0.0)
            fname = "DR_hydro_at_" + str(loc[i]) + ".csv"
            dataframe = np.transpose(np.array([np.array(time), np.array(time_sec), q, qc, qf, cc, cf]))
            np.savetxt(os.path.join(tdir,fname), dataframe, delimiter=",", \
                header="time, sec, q, qc, qf, cc, cf", fmt="%s, %i, %.5f, %.5f, %.5f, %.5f, %.5f")
            # png
            fig, ax = plt.subplots(3,1, figsize=(w_inchi,h_inchi))
            for ii in range(3):
                ax[ii].xaxis.set_major_formatter(DateFormatter('%b.%d\n%H:%M'))
            ax2=[ax[ii].twinx() for ii in range(3)]
            [ax2ax.xaxis.set_major_formatter(DateFormatter('%b.%d\n%H:%M')) for ax2ax in ax2]
            ax[0].plot(time,q,lw=0.5,color='k',ls='-',marker='',label='water & sed.')
            ax2[0].plot(time,qa.T[i],lw=0.5,color='k',ls='--',marker='',label='No.%d'%loc[i])
            ax[1].plot(time,qc,lw=0.5,color='k',ls='-',marker='',label='coarse sed.')
            ax2[1].plot(time,qca.T[i],lw=0.5,color='k',ls='--',marker='',label='No.%d'%loc[i])
            ax[2].plot(time,qf,lw=0.5,color='k',ls='-',marker='',label='fine sed.')
            ax2[2].plot(time,qfa.T[i],lw=0.5,color='k',ls='--',marker='',label='No.%d'%loc[i])
            ti=['(a)','(b)','(c)']
            for iii in range(3):
                # '3' means num of panels such as discharge, cc, cf
                # ax[i].set_xlim(datetime.datetime(1990,7,2,8), datetime.datetime(1990,7,2,14))
                ax[iii].set_xlim(time[0], time[-1])
                if iii==1: ax[iii].set_ylabel(u'Discharge (m$^3$/s) at {}'.format(str(loc[i])), fontsize=fsize)
                if iii<3: ax[iii].axes.xaxis.set_ticklabels([])
                ax[iii].set_ylim(0,)
                ax[iii].legend(loc='center right',fontsize=fsize)
                #ax2=ax[i].twinx()
                #ax2[i].yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
                ax2[iii].ticklabel_format(style='sci',axis='y',scilimits=(0,0))
                ax2[iii].set_ylim(0,)
                if i==0: ax2[iii].set_ylabel(u'Volume (m$^3$)',fontsize=fsize)
                ax2[iii].annotate('%s'%ti[i], xy=[0., 0.], xytext=[0.02,0.79], xycoords='axes fraction')
            #ax[0].axes.xaxis.set_ticklabels([])
            #ax[1].axes.xaxis.set_ticklabels([])
            #ax3=ax[2].twinx()
            #ax3.xaxis.set_major_formatter(DateFormatter('%b.%d\n%H:%M'))
            ax[2].set_xlabel(u'Date/Time {}'.format(time[0].year),fontsize=fsize)
            #plt.grid()
            plt.tight_layout()
            # fname = "DR_hydro_{}".format(str(loc[i])) + ".png"
            fname = os.path.splitext(fname)[0] + ".png"
            plt.savefig(os.path.join(tdir,fname),dpi=300)
            #plt.show()
            print("")
            print("make png file in {}".format(os.path.join(tdir, fname)))
            # https://qiita.com/Masahiro_T/items/037c118c9dd26450f8c8
            plt.clf()  # plt.savefig cause not to release memeory
            plt.close()
            # import psutil
            # mem = psutil.virtual_memory().free / 1e9
            # print(i, f'memory used: {mem} [GB]')

    # ti=['(a)','(b)','(c)']
    # for i in range(3):
    #     # ax[i].set_xlim(datetime.datetime(1990,7,2,8), datetime.datetime(1990,7,2,14))
    #     ax[i].set_xlim(time[0], time[-1])
    #     if i==1: ax[i].set_ylabel(u'Discharge (m$^3$/s) at {}'.format(str(i)), fontsize=fsize)
    #     if i<2: ax[i].axes.xaxis.set_ticklabels([])
    #     ax[i].set_ylim(0,)
    #     ax[i].legend(loc='center right',fontsize=fsize)
    #     #ax2=ax[i].twinx()
    #     #ax2[i].yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    #     ax2[i].ticklabel_format(style='sci',axis='y',scilimits=(0,0))
    #     ax2[i].set_ylim(0,)
    #     if i==1: ax2[i].set_ylabel(u'Volume (m$^3$)',fontsize=fsize)
    #     ax2[i].annotate('%s'%ti[i], xy=[0., 0.], xytext=[0.02,0.79], xycoords='axes fraction')
    # #ax[0].axes.xaxis.set_ticklabels([])
    # #ax[1].axes.xaxis.set_ticklabels([])
    # #ax3=ax[2].twinx()
    # #ax3.xaxis.set_major_formatter(DateFormatter('%b.%d\n%H:%M'))
    # ax[2].set_xlabel(u'Date/Time {}'.format(time[0].year),fontsize=fsize)
    # #plt.grid()
    # plt.tight_layout()
    # fname = "DR_hydro_{}".format(str(i)) + ".png"
    # plt.savefig(os.path.join(tdir,fname),dpi=300)
    # #plt.show()
    # print("")
    # print("make png file in {}".format(os.path.join(tdir, fname)))

import argparse
# parser.add_argument('--flag', action='store_true', required=True)
# parser.add_argument('--fruit', choices=['apple', 'banana', 'orange'])
# parser = argparse.ArgumentParser(\
    # prog="DR_out_hydro.py",
    # usage="python DR_out_hydro.py river_node_number_list",
    # description="description",
    # add_help=True,
    #)
parser = argparse.ArgumentParser()
tp = lambda x: list(map(int, x.split(",")))
parser.add_argument("nodes", \
                    type=tp, \
                    help="river nodes ids; ex) 1574 or 1574,2406", \
                    )
# parser.add_argument("-sdate", "--start_date", \
#                     type=str, \
#                     default="", \
#                     help="start date of x axis: format yyyy-m-d-H-M; ex. 1993-7-2-10-0", \
#                     )
# parser.add_argument("-edate", "--end_date", \
#                     type=str, \
#                     help="end date of x axis: format yyyy-m-d-H-M; ex. 1993-7-2-23-0", \
#                     )
parser.add_argument("-s", "--skip", \
                    type=int, \
                    default=1, \
                    help = "read files in 'output_DR' by each SKIP interval"
                    )
parser.add_argument("--version", version="%(prog)s 1.0.0",
                    action="version",
                    default=False)
args = parser.parse_args()

if len(args.nodes) > 0 and args.skip > 0:
    make_DR_fig(loc=args.nodes, SKIP=args.skip)
else:
    # https://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu
    import sys
    parser.print_help(sys.stderr)
    sys.exit(1)

