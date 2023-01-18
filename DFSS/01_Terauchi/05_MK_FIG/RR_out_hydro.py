import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import datetime
from matplotlib.dates import DateFormatter
import sys

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
fname='Terauchi_dam.csv'
dt=3600
data=np.loadtxt(fname,dtype=str,delimiter=',',skiprows=0)
iskp=0
q_obs=np.transpose(data)[1].astype('float')
cont=np.transpose(data)[0][0]
date,time=cont.split()[0].split('/'),cont.split()[1].split(':')
year,month,day=int(date[0]),int(date[1]),int(date[2])
hour,minute=int(time[0]),int(time[1])
time0=datetime.datetime(year,month,day,hour,minute)
time0=time0-datetime.timedelta(seconds=0)
time=[time0+datetime.timedelta(seconds=dt*i) for i in range(len(q_obs))]
plt.plot(time,q_obs,lw=0.6,ls='--',color='k',            label=r'obs.')
#################################################################################

ls='\n'
dirs=glob.glob(os.path.join('..','04_EXE_SIMU','01_RR_DR'))
#loc=[4575,4741] # 343:St.1, 1574 1562:st.7
loc=sys.argv[1:]
loc=[int(s) for s in loc]
print('out @',loc)
model='RR'
lgs=[[1] for i in range(len(dirs))]
for lg,tdir in zip(lgs,dirs):
    if model=='RR': files=glob.glob(os.path.join(tdir,'output_RR','res_1d_*.txt'))
    files.sort()
    iend=len(files)
    label=tdir.split('_')[-1]
    fn=os.path.join(tdir,'output_misc','basin_rain.txt')
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
    for i0,afile in enumerate(files[:]):
        res=np.loadtxt(afile,usecols=(3,4,5,6))
        sec=float(os.path.basename(afile).split('_')[2][:-4])
        #print (afile)
        z=np.transpose(res)[0]
        h=np.transpose(res)[1]
        u=np.transpose(res)[2]
        q=np.transpose(res)[3]
        hydro.append([q[l-1] for l in loc])
        hhh.append([h[l-1] for l in loc])
        uuu.append([u[l-1] for l in loc])
        time.append(time0+datetime.timedelta(seconds=sec))
        #print(sec,time[-1])
    for i in range(len(loc)):
        plt.plot(time,np.transpose(hydro)[i],lw=lg[0],marker='',label='sim.')
plt.xlim(datetime.datetime(2017,7,5,0), datetime.datetime(2017,7,7,0))
plt.xlabel(u'Year (2017)',fontsize=fsize)
plt.ylabel(u'Discharge (m$^3$/s)',fontsize=fsize)
plt.ylim(0,)
plt.legend(fontsize=fsize)
#plt.grid()
plt.tight_layout()
plt.savefig('%s_hydro.png'%model,dpi=300)
#plt.show()




