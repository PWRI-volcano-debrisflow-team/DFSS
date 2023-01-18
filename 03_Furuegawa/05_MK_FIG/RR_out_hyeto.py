import os
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
from matplotlib.font_manager import FontProperties
fp = FontProperties(fname=r'C:\WINDOWS\Fonts\msgothic.ttc')
##################################
fig, ax = plt.subplots(1, 1, figsize=(w_inchi,h_inchi))
ax.xaxis.set_major_formatter(DateFormatter('%b%d\n%H:%M\n'))
#ax.xaxis.set_major_formatter(DateFormatter('%H'))
ax2=ax.twinx()
ax2.xaxis.set_major_formatter(DateFormatter('%b.%d\n%H:%M'))
ls='\n'
fs=','
tdir=os.path.join('..','04_EXE_SIMU','01_RR_DR')
fn=os.path.join(tdir,'output_misc','basin_rain.txt')
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
ax.plot(datetime.datetime(1990,7,2,9,20),80,marker='^',color='k')
arain=[sum(hyet[0:i]) for i in range(1,len(hyet)+1)]
ax2.plot(time1,arain,lw=0.8,ls='-',color='k')
ax2.set_ylim(0,500)
#plt.xlim(datetime.datetime(2017,7,5,0), datetime.datetime(2017,7,6,12))
ax.set_xlabel(u'Year (1990)',fontsize=fsize)
plt.xticks(fontsize=fsize*0.8)
ax.set_ylim(100,0)
#ax.set_ylabel(u'%d分雨量 (mm)'%int(dt.seconds/60),fontproperties=fp,fontsize=fsize)
#ax2.set_ylabel(u'累積雨量 (mm)',fontproperties=fp,fontsize=fsize)
ax.set_ylabel(u'%d mih. rainfall (mm)'%int(dt.seconds/60),fontsize=fsize)
ax2.set_ylabel(u'Accumulated rainfall (mm)',fontsize=fsize)
#plt.legend(fontsize=fsize)
plt.tight_layout()
#plt.grid()
#plt.yscale("log")
plt.savefig('RR_hyeto.png',dpi=300)
#plt.show()
