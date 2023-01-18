import numpy as np

qnum=6
qdt=600
flux=0.01 # flox at inlet
ista=274 # loc on inlet

tend=qnum*qdt
ista=ista-1
ls='\n'
fname='streamConfiguration_inRR.txt'
data=np.loadtxt(fname,skiprows=2,dtype=str,delimiter=',')
i_2d=data.T[14].astype(int)
j_2d=data.T[15].astype(int)
itr=int(tend/qdt)
ts=np.linspace(0,tend,itr+1)
data=[0 for i in range(len(i_2d))]
f=open('flowVolume_RRtoDR.txt','w')
for ij in [i_2d,j_2d]:
    s='{:>20s}'.format('i_2d')+''.join(['{:12d}'.format(i) for i in ij[:]])
    f.write(s+ls)
for t in ts:
    v=t*flux
    data[ista]=v
    s='  time(s)={:>10d}'.format(int(t))+''.join(['{:12.3e}'.format(f) for f in data[:]])
    f.write(s+ls)
f.close()
