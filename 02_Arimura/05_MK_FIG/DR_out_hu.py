import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import sys
#####################################################
r_gold=(1.+5.**0.5)*0.5
r_silver=(1.+2.**0.5)
r_yamato=2.**0.5
w_inchi=3.14961 # 80mm
h_inchi=w_inchi/r_gold
r_top=0
r_bottom=0
r_left=0
r_right=0
fsize=8
#####################################################
plt.figure(figsize=(w_inchi,h_inchi))
plt.rcParams["font.size"] = fsize
#####################################################
fname='201408291557.txt'
iskp=-8
data=np.loadtxt(fname,dtype=str,delimiter=',',skiprows=1)
h_obs=np.transpose(data)[2].astype('float')
v_obs=np.transpose(data)[1].astype('float')
time=[i*1+iskp for i in range(len(h_obs))]
v_obs=[f for f in v_obs]
plt.plot(time,h_obs,lw=0.6,ls='--',color='k',            label=r'h (obs.)')
plt.plot(time[:-1],v_obs[1:],lw=0.6,ls='--',color='gray',label=r'u (obs.)')
#####################################################
rm_obs=np.transpose(data)[3].astype('float')
import math
tant=math.tan(math.radians(4))
tanp=math.tan(math.radians(37))
sig=2.65
g=9.81
#def c_rho(rhm,tant,tanp):
#    return rhm*(tanp-tant)/tanp
#def c_cc(sig,rho,tant,tanp):
#    return tant/(tanp-tant)/(sig/rho-1.)
#rho=[c_rho(x1/g,tant,tanp) for x1 in rm_obs]
#cc=[c_cc(sig,x1,tant,tanp)for x1 in rho]
#print cc
#plt.plot(time[:-1],cc[1:],'r--',label='c(obs.)')
######################################################
loc=sys.argv[1:]
loc=[int(s) for s in loc]
print('out @',loc)
dirs=[os.path.join('..','04_EXE_SIMU','01_RR_DR')]
for tdir in dirs:
    files=glob.glob(os.path.join(tdir,'output_DR','res_??????????.txt'))
    files.sort()
    for l in loc:
        chk=[np.loadtxt(afile,skiprows=1)[l-1] for afile in files]
        # ---> set time step to obs (1 min)
        chk=[sum(chk[i:i+6])/6. for i in range(0,len(chk),6)]
        time=[i for i in range(len(chk))]
        # <--- set time step to obs (1 min)
        chk=np.array(chk).T
        x=[np.array(chk[i]) for i in range(9)]
        # z,h,u,q,cc,cf,rho,e,z0
        # 1 2 3 4  5  6   7 8  9
        h=x[2]
        u=x[3]
        q=x[4]
        plt.plot(time,h,lw=0.6,ls='-',color='k',label=r'h (sim.)')
        plt.plot(time,u,lw=0.6,ls='-',color='gray',label=r'u (sim.)')
plt.xlim(0,25)
plt.ylim(0,)
plt.yticks([1,2,3,4,5])
plt.xlabel('Time (min.)')
plt.ylabel(r'$\it{h}$ (m), $\it{u}$ (m/s)')
#plt.yscale('log')
plt.tight_layout()
plt.legend(loc='upper right',fontsize=fsize*0.8)
plt.savefig('DR_hu.png',dpi=300)
#plt.show()
