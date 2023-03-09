def make_RR_fig(EVENT, RAIN=80):
    import os, glob, shutil
    import numpy as np
    import matplotlib.pyplot as plt
    import datetime
    from matplotlib.dates import DateFormatter
    ##################################
    r_gold=(1.+5.**0.5)*0.5
    r_silver=(1.+2.**0.5)
    r_yamato=2.**0.5
    w_inchi=4
    h_inchi=w_inchi/r_gold*0.8
    r_top=0
    r_bottom=0
    r_left=0
    r_right=0
    fsize=7.5
    ##################################
    #plt.figure(figsize=(w_inchi,h_inchi))
    plt.rcParams["font.size"] = fsize
    # japanese font setting
    # https://qiita.com/yniji/items/3fac25c2ffa316990d0c
    # from matplotlib.font_manager import FontProperties
    # fp = FontProperties(fname=r'C:\WINDOWS\Fonts\msgothic.ttc')
    ##################################
    tdir = "output_RR_fig"
    file_list = glob.glob("*")
    if tdir in file_list:
        tdir = os.path.join(os.getcwd(), tdir)
    else:
        tdir = os.path.join(os.getcwd(), tdir)
        os.mkdir(tdir)
    ##################################
    fig, ax = plt.subplots(1, 1, figsize=(w_inchi,h_inchi))
    ax.xaxis.set_major_formatter(DateFormatter('%b%d\n%H:%M\n'))
    #ax.xaxis.set_major_formatter(DateFormatter('%H'))
    ax2=ax.twinx()
    ax2.xaxis.set_major_formatter(DateFormatter('%b.%d\n%H:%M'))
    ls='\n'
    fs=','
    #srcdir=os.path.join('..','04_EXE_SIMU','01_RR_DR')
    srcdir = os.getcwd()
    fn=os.path.join(srcdir,'output_RR_misc','basin_rain.txt')
    shutil.copy(fn, os.path.join(tdir, "basin_rain.txt"))
    txt=np.loadtxt(fn,dtype=str,delimiter=',',skiprows=1)
    date=txt.T[0]
    hyet=txt.T[1].astype(float)
    time1=[]
    for line in txt.T[0]:
        date=line.split()[0].split('/')
        year,month,day=int(date[0]),int(date[1]),int(date[2])
        time=line.split()[1].split(':')
        hour,minute=int(time[0]),int(time[1])
        time0=datetime.datetime(year,month,day,hour,minute)
        time1.append(time0)
    dt=time1[1]-time1[0]
    ax.fill_between(time1,hyet,linewidth=0.,facecolor='gray',alpha=1.,step='pre')
    #ax.plot(datetime.datetime(1990,7,2,9,20),80,marker='^',color='k')
    # draw triangle of debris flow events
    #ax.plot(time1[0],80,marker='^',color='k')
    if "noevent" in EVENT:
        # user input no event
        pass
    else:
        timing = []
        for event in EVENT:
            # user input a event
            year, month, day, hour, min = [int(e) for e in event.split("-")]
            timing.append(datetime.datetime(year, month, day, hour, min))
        rain_intensity = [RAIN for i in range(len(timing))]
        print(timing, rain_intensity)
        # in plot_date, marker is designated as fmt keyword.
        ax.plot_date(timing,rain_intensity,fmt='^',color='k')
    #
    arain=[sum(hyet[0:i]) for i in range(1,len(hyet)+1)]
    ax2.plot(time1,arain,lw=0.8,ls='-',color='k')
    ax2.set_ylim(0,500)
    #plt.xlim(datetime.datetime(2017,7,5,0), datetime.datetime(2017,7,6,12))
    yearlabel = "{}".format(time1[0].year)
    ax.set_xlabel(yearlabel,fontsize=fsize)
    plt.xticks(fontsize=fsize*0.8)
    ax.set_ylim(100,0)
    #ax.set_ylabel(u'%d分雨量 (mm)'%int(dt.seconds/60),fontproperties=fp,fontsize=fsize)
    #ax2.set_ylabel(u'累積雨量 (mm)',fontproperties=fp,fontsize=fsize)
    ax.set_ylabel(u'%d min. rainfall (mm)'%int(dt.seconds/60),fontsize=fsize)
    ax2.set_ylabel(u'Accumulated rainfall (mm)',fontsize=fsize)
    #plt.legend(fontsize=fsize)
    plt.tight_layout()
    #plt.grid()
    #plt.yscale("log")
    plt.savefig(os.path.join(tdir,'RR_hyeto.png'),dpi=300)
    #plt.show()


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-e", "--event", \
                    type=str, \
                    default="noevent", \
                    help = """draw the debris-flow timing as a triangle mark; format: yyyy-m-d-H-M;
                        Ex1 one event) 1990-7-2-9-20 OR Ex2 twe events) 1990-7-2-9-20,1990-7-2-10-20
                        """
                    )
parser.add_argument("--version", version="%(prog)s 1.0.0",
                    action="version",
                    default=False)
args = parser.parse_args()

event = args.event.split(",")
if len(event) > 0:
    make_RR_fig(EVENT=event)
else:
    import sys
    parser.print_help(sys.stderr)
    sys.exit(1)


