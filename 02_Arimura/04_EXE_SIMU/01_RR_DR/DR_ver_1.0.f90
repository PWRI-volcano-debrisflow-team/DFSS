! Program: DR_ver_1.0.f90
! Copyright (2021) by Public Works Research Institute (P.W.R.I.)
! License: CC-BY-SA
!-------------------------------------------------------------------------------
!        from  i    to
! FX      q    q    q
!        ->   ->   ->
! CV | z  | z  | z  |
!    |from| i  | to |
!-------------------------------------------------------------------------------
implicit none
real(8), parameter :: PI=acos(-1.d0),deg2rad=PI/180.d0,rad2deg=180.d0/PI, &
                      g=9.8d0,dtlim=0.001d0,cfl=0.01d0
real(8), allocatable :: zini(:),z(:),zs(:),b(:),h(:),ho(:),um(:), &
                        cc(:),cco(:),cf(:),cfo(:),e(:),rho(:),rhm(:),tant_e(:),n(:), &
                        q_in(:),q_out(:),qcc_in(:),qcc_out(:),qcf_in(:),qcf_out(:), &
                        hs(:),d(:),c(:),pw(:),sf(:),hsc(:)
real(8), allocatable :: u(:),uo(:),grad_z(:),grad_w(:),grad_s(:),q(:),qcc(:),qcf(:), &
                        qs(:),qs_in(:),qs_out(:),q_a(:),qc_a(:),qf_a(:)
integer :: ni_inlet,ni_outlet
integer, allocatable :: i_inlet(:),i_outlet(:)
real(8), allocatable :: qin0(:),cin0(:),rin0(:),fsr_in0(:,:),fsr_in(:)
real(8) :: qin,cin,rin
real(8) :: dt,time,end_time,dout_time,fout_time,qdt
integer :: i,iend,k,itsp,itsp_end,itsp_dout,itsp_fout,qnum,itsp_q,i_umax,itmp
integer :: ni_scal,ni_vect
integer, allocatable :: i_node(:),i_to(:),i_from(:,:),i_scal(:),i_vect(:),ifld(:)
real(8), allocatable :: dx_scal(:),dx_vect(:)
integer :: iflg_1d
real(8) :: ba,d_1d,rn_1d,wid_m,wid_p,wid_min,dep_m,dep_p,dep_min
real(8) :: sig,dm,tanp,csta,cstad,pc,pf,rhow,gam,cmin,hmin,hlim,umax,humax, &
           ks,lam,fs2,cohe,hs_ini
real(8) :: cin_sum,cout_sum,dz_sum,q_sum,qcc_sum,qcf_sum,qs_sum, &
           h_sum,hs_sum,cc_sum,cf_sum,fsr_sum
logical, allocatable :: mask(:)
character(len=100) :: fname,out_dir,ctmp
integer :: itr
real(8) :: flmax,time0,time1,dxmin,dt0,dtmin
! out ---> 2D
real(8), allocatable :: z_2d(:),va_out(:),vc_out(:),vf_out(:)
integer, allocatable :: iacc(:)
integer :: iaccmax,iacc_up
!
open(1,file='DR_input.txt')
read(1,*)
read(1,*)iend
read(1,*)iacc_up
read(1,*)
read(1,*)dt
read(1,*)dout_time
read(1,*)fout_time
read(1,*)
read(1,*)qnum
read(1,*)qdt
!
allocate(zini(iend),z(iend),zs(iend),b(iend),h(iend),ho(iend),um(iend), &
         cc(iend),cco(iend),cf(iend),cfo(iend),e(iend),rho(iend),rhm(iend),tant_e(iend),n(iend), &
         q_in(iend),q_out(iend),qcc_in(iend),qcc_out(iend),qcf_in(iend),qcf_out(iend), &
         hs(iend),d(iend),c(iend),pw(iend),sf(iend),hsc(iend))
allocate(u(iend),uo(iend),grad_z(iend),grad_w(iend),grad_s(iend),q(0:iend),qcc(0:iend),qcf(0:iend), &
         qs(0:iend),qs_in(iend),qs_out(iend),q_a(0:iend),qc_a(0:iend),qf_a(0:iend))
allocate(qin0(0:qnum),cin0(0:qnum),rin0(0:qnum),fsr_in0(0:qnum,iend),fsr_in(iend))
allocate(i_node(iend),i_to(iend),i_from(iend,3),dx_scal(iend),dx_vect(iend))
allocate(mask(iend),ifld(iend),iacc(iend),z_2d(iend),va_out(iend),vc_out(iend),vf_out(iend))
!
!open(10,file='bound_BSR.txt')
do i=0,qnum-1 ! i=0: initial
!  read(10,*)qin0(i),cin0(i),rin0(i)
  qin0(i)=0.d0
  cin0(i)=0.d0
  rin0(i)=0.d0
  rin0(i)=rin0(i)/dble(1000*3600)
end do
qin0(qnum)=qin0(qnum-1)
cin0(qnum)=cin0(qnum-1)
rin0(qnum)=rin0(qnum-1)
!close(10)
read(1,*)
read(1,*)dm
read(1,*)pf
read(1,*)ks ! (mm/h)
ks=ks/100.0d0 ! (mm/h) -> (m/s)
read(1,*)csta
read(1,*)rhow
read(1,*)sig
read(1,*)tanp ! deg
tanp=tan(tanp*deg2rad) ! deg ---> rad
read(1,*)cohe
read(1,*)iflg_1d,hs_ini,fname
if(iflg_1d/=0)then
  open(10,file=fname)
  read(10,*)
  do i=1,iend
    read(10,*)itmp,hs(i)
  end do
  close(10)
else
  hs(:)=hs_ini
end if
read(1,*)
read(1,*)hlim
read(1,*)
read(1,*)iflg_1d
read(1,*)fname
!
open(10,file=fname)
read(10,*)i
if(i/=iend)then
write(*,*)'Check iend'
end if
read(10,*)
do i=1,iend
  read(10,*)i_node(i),i_to(i),i_from(i,1),i_from(i,2),i_from(i,3),&
  z_2d(i),z(i),zs(i),dx_scal(i),b(i),n(i),iacc(i)
  dx_vect(i)=dx_scal(i)
end do
close(10)
dxmin=minval(dx_vect(:))
!
read(1,*)d_1d
read(1,*)rn_1d
read(1,*)wid_m
read(1,*)wid_p
read(1,*)wid_min
read(1,*)dep_m
read(1,*)dep_p
read(1,*)dep_min
close(1)
!
if(iflg_1d==0)then
  write(*,*)'Stream config. <--- DR_input.txt'
  do i=1,iend
    ba=dxmin*dxmin*float(iacc(i))/1000000.d0
    b(i)=max(wid_min,wid_m*ba**wid_p)
    d(i)=max(dep_min,dep_m*ba**dep_p)
    z(i)=z_2d(i)-d(i)
    n(i)=rn_1d
    zs(i)=z(i)-d_1d
  end do
else
  write(*,*)'Stream config. <---',trim(adjustl(fname))
end if
zini(:)=z(:)
ifld(:)=0
!
open(20,file='DR_streamChecker.txt')
write(20,'(a5,",",*(a10,:,","))')'id','z_2d','z_1d','zs_1d','dx','b','n'
do i=1,iend
  write(20,'(i5,",",5(f10.3,:,","),*(f10.5,:,","))')i_node(i),z_2d(i),z(i),zs(i),dx_scal(i),b(i),n(i)
end do
close(20)
!
open(10,file='flowVolume_RRtoDR.txt')
read(10,*)
read(10,*)
do i=0,qnum-1 ! i=0: initial
  read(10,'(a10,i10,*(e12.2))')ctmp,k,fsr_in0(i,:)
  fsr_in0(i,:)=fsr_in0(i,:)/(b(:)*dx_scal(:))/qdt
end do
fsr_in0(1:qnum-1,:)=fsr_in0(1:qnum-1,:)-fsr_in0(0:qnum-2,:)
fsr_in0(qnum,:)=fsr_in0(qnum-1,:)
where(fsr_in0<0.d0)fsr_in0=0.d0
close(10)
!do i=0,qnum-1 ! i=0: initial
!  write(*,'(f5.1,a,f10.3,a)')dble(i)*qdt/3600.d0,'(h)',&
!  sum(fsr_in0(i,:))/dble(size(fsr_in0(i,:)))*(1000.d0*3600.d0),' (mm/h)'
!end do
!
mask=i_from(:,1)/=0
ni_scal=count(mask)
allocate(i_scal(ni_scal))
i_scal(1:ni_scal)=pack(i_node(:),mask)
!
mask=i_to(:)/=0
ni_vect=count(mask)
allocate(i_vect(ni_vect))
i_vect(1:ni_vect)=pack(i_node(:),mask)
!
mask=i_from(:,1)==0
ni_inlet=count(mask)
allocate(i_inlet(ni_inlet))
i_inlet(1:ni_inlet)=pack(i_node(:),mask)
!
mask=i_to(:)==0
ni_outlet=count(mask)
allocate(i_outlet(ni_outlet))
i_outlet(1:ni_outlet)=pack(i_node(:),mask)
!
out_dir='output_DR/'
write(*,*)trim(adjustl(out_dir))

lam=1.d0-csta
fs2=0.d0
iaccmax=iacc(iacc_up)
c(:)=cohe
cstad=csta
pc=1.d0-pf
hmin=0.01d0
cmin=0.001d0
gam=1.d0
!
end_time=dble(qnum-1)*qdt
itsp_end=nint(end_time/dt)
itsp_dout=nint(dout_time/dt)
itsp_fout=nint(fout_time/dt)
itsp_q=nint(qdt/dt)
!
h(:)=0.d0
u(:)=0.d0
e(:)=0.d0
ho(:)=h(:)
uo(:)=u(:)
cc(:)=0.d0
cf(:)=0.d0
cco(:)=cc(:)
cfo(:)=cf(:)
q(:)=0.d0
qcc(:)=0.d0
qcf(:)=0.d0
dz_sum=0.d0 !sum(b(1:iend-1)*dx*(z(1:iend-1)-zs(1:iend-1))*csta)
cin_sum=0.d0
cout_sum=0.d0
!hs(:)=0.d0
pw(:)=0.d0 ! <--- dummy
d(:)=z(:)-zs(:)
!
q_sum=0.d0
qcc_sum=0.d0
qcf_sum=0.d0
qs_sum=0.d0
h_sum=0.d0
hs_sum=0.d0
cc_sum=0.d0
cf_sum=0.d0
fsr_sum=0.d0
!
q_a(:)=0.d0
qc_a(:)=0.d0
qf_a(:)=0.d0
!
umax=0.d0
humax=0.d0
i_umax=0
!
va_out(:)=0.d0
vc_out(:)=0.d0
vf_out(:)=0.d0
!
! ---> cal start
open(200,file='DR_summary.txt')
write(200,'(a10,*(a15,:))')'time(s)','q_sum','qs_sum','h_sum','hs_sum', &
                           'qcc_sum','qcf_sum','cc_sum','cf_sum','dz_sum','fsr_sum'
open(201,file='flowVolume_ws_DRtoDF.txt')
open(202,file='flowVolume_cs_DRtoDF.txt')
open(203,file='flowVolume_fs_DRtoDF.txt')
dt0=dt
dtmin=dt
time=0.d0
time0=0.d0
write(201,'(a10,"",i10,*(E20.10e3,:,""))')'time(s)=',idnint(time),va_out(:)
write(202,'(a10,"",i10,*(E20.10e3,:,""))')'time(s)=',idnint(time),vc_out(:)
write(203,'(a10,"",i10,*(E20.10e3,:,""))')'time(s)=',idnint(time),vf_out(:)
!
do itsp=1,itsp_end
  time=time0+dt
  k=itsp/itsp_q
  qin=qin0(k)+(qin0(k+1)-qin0(k))*(itsp-k*itsp_q)/dble(itsp_q)
  cin=cin0(k)+(cin0(k+1)-cin0(k))*(itsp-k*itsp_q)/dble(itsp_q)
  fsr_in(:)=fsr_in0(k,:)+(fsr_in0(k+1,:)-fsr_in0(k,:))*(itsp-k*itsp_q)/dble(itsp_q)
  rin=rin0(k+1)
  !
  flmax=0.d0
  call cal_u(ni_vect,i_vect,iend,3,i_from,i_to,dtlim,g,n,z,b,h,ho,grad_w,dx_vect,&
       u,uo,hlim,rho,rhm,sig,dm,csta,tanp,cc,umax,i_umax,humax)
  flmax=umax
  !
  time1=time0
  itr=1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ((dxmin*cfl) < flmax)then
      itr=max(1,ceiling(flmax/(dxmin*cfl)))
      dt=max(dtlim,dt0/float(itr))
      !write(*,'(a,f8.4)', advance='yes')'dt=',dt
      if(dt < dtmin)then
        dtmin=dt
        itr=nint(dt0/dt)
      end if
    end if
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  umax=0.d0
  i_umax=0
  humax=0.d0
  do while (itr>0)
    d(:)=z(:)-zs(:)
    call cal_hs(iend,i_node,iend,dt,dx_scal,fsr_in,lam,b,hs,iend,qs_in,qs_out)
    where(hs(:)>d(:))
      h(:)=h(:)+(hs(:)-d(:))*lam
      hs(:)=d(:)
    end where
    !
    call cal_sf(ni_scal,i_scal,iend,g,rhow,sig,tanp,lam,hlim,d,grad_z,hs,h,pw,c,sf,hsc)
    where(1.d0<sf(:) .and. d(:)==hs(:))
      ifld(:)=1
    end where
    do i=1,ni_inlet
      k=i_to(i_inlet(i))
      do while (k/=0)
        if( ifld(i_from(k,1))==1 )then
          ifld(k)=1
        end if
        k=i_to(k)
      end do
    end do
    where(ifld==0)e=0.d0
    !
    call cal_h(iend,i_node,iend,dt,dx_scal,b,h,iend,q_in,q_out,e,csta,cstad)
    call cal_cc(iend,i_node,iend,dt,dx_scal,b,h,ho,cc,cco,iend,qcc_in,qcc_out,e,gam,pc)
    call cal_cf(iend,i_node,iend,dt,dx_scal,b,h,ho,cf,cfo,cc,cco,iend,qcf_in,qcf_out,e,gam,pf,cstad)
    call cal_z(iend,i_node,iend,dt,csta,z,zs,e,tant_e)
    z(i_outlet(1))=zs(i_outlet(1))+(z(i_from(i_outlet(1),1))-zs(i_from(i_outlet(1),1)))
    !
    call cal_vout(iend,i_node,iend,dt,dx_scal,b,h,z,z_2d,cc,cf,iacc,iaccmax, &
                  va_out,vc_out,vf_out)
    !
    call cal_qs(ni_vect,i_vect,iend,i_to,b,hs,grad_s,qs,ks)
    i=i_from(i_outlet(1),1)
    qs(i_outlet(1))=b(i)*ks*hs(i)*min(1.0d0,abs(grad_s(i)))
    qs_in(i_scal(:))=qs(i_from(i_scal(:),1))+qs(i_from(i_scal(:),2))+qs(i_from(i_scal(:),3))
    qs_out(:)=qs(i_node(:))
    do i=1,ni_inlet
      qs_in(i_inlet(i))=0.d0
    end do
    !
    call cal_u(ni_vect,i_vect,iend,3,i_from,i_to,dt,g,n,z,b,h,ho,grad_w,dx_vect,&
         u,uo,hlim,rho,rhm,sig,dm,csta,tanp,cc,umax,i_umax,humax)
    u(i_outlet(1))=u(i_from(i_outlet(1),1))
    !
    call cal_q(ni_vect,i_vect,iend,i_to,b,h,u,q(0:iend),cc,cf,qcc(0:iend),qcf(0:iend))
    ! one outlet
    i=i_outlet(1)
    q(i)=b(i)*h(i)*u(i)
    qcc(i)=q(i)*cc(i)
    qcf(i)=q(i)*cf(i)*(1.d0-cc(i))
    !
    q_in(i_scal(:))=q(i_from(i_scal(:),1))+q(i_from(i_scal(:),2))+q(i_from(i_scal(:),3))
    qcc_in(i_scal(:))=qcc(i_from(i_scal(:),1))+qcc(i_from(i_scal(:),2))+qcc(i_from(i_scal(:),3))
    qcf_in(i_scal(:))=qcf(i_from(i_scal(:),1))+qcf(i_from(i_scal(:),2))+qcf(i_from(i_scal(:),3))
    !
    q_out(:)=q(i_node(:))
    qcc_out(:)=qcc(i_node(:))
    qcf_out(:)=qcf(i_node(:))
    !
    do i=1,ni_inlet
      q_in(i_inlet(i))=qin
      qcc_in(i_inlet(i))=qin*pc*cin
      qcf_in(i_inlet(i))=qin*pf*cin
    end do
    !
    rho(i_node)=(sig-rhow)*cf(:)+rhow
    rhm(:)=(sig-rho(:))*cc(:)+rho(:)
    !
    do k=1,ni_vect
      i=i_vect(k)
      grad_w(i)=(z(i_to(i))+h(i_to(i))-z(i)-h(i))/dx_scal(i)
      grad_s(i)=(zs(i_to(i))+hs(i_to(i))-zs(i)-hs(i))/dx_scal(i)
      grad_z(i)=(z(i_to(i))-z(i))/dx_scal(i)
    end do
    grad_w(i_outlet(1))=grad_w(i_from(i_outlet(1),1))
    grad_s(i_outlet(1))=grad_s(i_from(i_outlet(1),1))
    grad_z(i_outlet(1))=grad_z(i_from(i_outlet(1),1))
    !
    um(i_scal(:))=(u(i_scal(:))+u(i_from(i_scal(:),1)))*0.5d0
    where(h(i_inlet(:))<hlim)
      um(i_inlet(:))=u(i_inlet(:))*0.5d0
    elsewhere
      um(i_inlet(:))=(u(i_inlet(:))+qin/(h(i_inlet(:))*b(i_inlet(:))))*0.5d0
    end where
    um(i_inlet(:))=u(i_inlet(:))
    !
    call cal_e(iend,i_node,iend,dt,dx_scal,b,z,zs,sig,rho,csta,cstad,pf,pc,tanp,cc,cf,e,h,um,grad_z, &
               q_in,q_out,qcc_in,qcc_out,qcf_in,qcf_out,tant_e)
    !
    uo(:)=u(:)
    ho(:)=h(:)
    cco(:)=cc(:)
    cfo(:)=cf(:)
    !
    q_a(:)=q_a(:)+q(:)*dt
    qc_a(:)=qc_a(:)+qcc(:)*dt
    qf_a(:)=qf_a(:)+qcf(:)*dt
    !
    cin_sum=cin_sum+sum(qcc_in(i_inlet(:)))*dt+sum(qcf_in(i_inlet(:)))*dt
    cout_sum=cout_sum+qcc(i_outlet(1))*dt+qcf(i_outlet(1))*dt
    q_sum=q_sum+q(i_outlet(1))*dt
    qcc_sum=qcc_sum+qcc(i_outlet(1))*dt
    qcf_sum=qcf_sum+qcf(i_outlet(1))*dt
    qs_sum=qs_sum+qs(i_outlet(1))*dt
    h_sum=sum(h(:)*b(:)*dx_scal(:))
    hs_sum=sum(hs(:)*b(:)*dx_scal(:))*lam
    dz_sum=sum(d(:)*b(:)*dx_scal(:))
    cc_sum=sum(h(:)*cc(:)*b(:)*dx_scal(:))
    cf_sum=sum(h(:)*(1.d0-cc(:))*cf(:)*b(:)*dx_scal(:))
    fsr_sum=fsr_sum+sum(fsr_in(:))*dt
    !
    time1=time1+dt
    itr=itr-1
  end do ! do while
  !
  time0=time
  dt=dt0
  ! output
  if (mod(itsp, itsp_dout) ==0) then
    write(*,'(a,i8,5(a,f10.3))')'time(s)=',nint(time), &
    ', qin=',sum(q_in(i_inlet(:))),', qout=',q(i_outlet(1)),',   umax=',umax
    write(*,'(a,f8.3,2(a,f10.3),(a,i6))')'time(h)=',(time/3600.d0), &
    '  qcc=',sum(qcc_in(i_inlet(:))),', qout=',qcc(i_outlet(1)),', i@umax=',i_umax
    write(*,'(a,5(a,f10.3))')'                ', &
    '  qcf=',sum(qcf_in(i_inlet(:))),', qout=',qcf(i_outlet(1)),', h@umax=',humax
    write(*,'(2(a,f10.5))')'dt_min=',dtmin,', dt_lim=',dtlim
    umax=0.d0
    humax=0.d0
    dtmin=dt
    ! ---> summary
    write(200,'(i10,"",*(e15.5e3,:,""))')nint(time),q_sum,qs_sum,h_sum,hs_sum, &
                                       qcc_sum,qcf_sum,cc_sum,cf_sum,dz_sum,fsr_sum
    ! <--- summary
    write(201,'(a10,"",i10,*(E20.10e3,:,""))')'time(s)=',idnint(time),va_out(:)
    write(202,'(a10,"",i10,*(E20.10e3,:,""))')'time(s)=',idnint(time),vc_out(:)
    write(203,'(a10,"",i10,*(E20.10e3,:,""))')'time(s)=',idnint(time),vf_out(:)
  end if
  ! ---> output to files
  if(mod(itsp, itsp_fout) ==0) then
    write(fname,'(i10.10)')idnint(time)
    fname=trim(adjustl(out_dir))//'res_'//trim(adjustl(fname))//'.txt'
    open(100,file=fname)
    write(100,'(a5,15a12,3a15)')'id','z','h','u','q','cc','cf','rho','e',&
    'zini','grd_w','hs','qs','sf','hsc/d','ifld','qa','qcca','qcfa'
    do i=1,iend
      write(100,'(i5,14f12.5,i12,3e15.5e3)')i,z(i),h(i),um(i),q(i),cc(i),cf(i),rho(i),e(i),zini(i),&
      grad_w(i),hs(i),qs(i),sf(i),hsc(i),ifld(i),q_a(i),qc_a(i),qf_a(i)
    end do
    close(100)
  end if
  ! <--- output to files
end do
! <--- cal start
close(200)
close(201)
close(202)
close(203)
write(*,*)'--- nomal end ---'
!
stop
end
!
!
!
subroutine cal_u(ni_vect,i_vect,iend,nifrom,i_from,i_to, dt,g,n,z,b,h,ho,grad,&
                 dx,u,uo,hlim,rho,rhm,sig,dm,csta,tanp,c,umax,i_umax,humax)
use omp_lib
implicit none
real(8), parameter :: PI=acos(-1.d0),deg2rad=PI/180.d0,rad2deg=180.d0/PI, &
                      hm=5.d0/3.d0,rm=4.d0/3.d0
real(8), parameter :: f1d5=1.d0/5.d0,f1d3=1.d0/3.d0,f5d3=5.d0/3.d0,f2d3=2.d0/3.d0, &
                      kf=0.16d0,kd=0.0828d0,e2=0.85d0**2.d0
integer :: i,iend,j,i_umax,iloc(1),ni_vect,i_vect(ni_vect),nifrom,i_from(iend,nifrom),i_to(iend)
real(8) :: dt,g,x1,x2,x3,hlim,h(iend),ho(iend),u(iend),uo(iend), &
           h1,b1,c1,rho1,rhm1,n1,z1,h2,z2,hnw, &
           c(iend),rho(iend),rhm(iend),z(iend),b(iend),n(iend), &
           sig,csta,dm,taub,tauy,fb,fb1,fb2,tanp,fd,ff,cost, &
           grad(iend),dx(iend),umax,uhb(nifrom),humax,htmp(iend)
!
umax=0.d0
u(:)=0.d0
htmp(:)=0.d0
humax=0.d0
i_umax=0
!$omp parallel do private(i,j,cost,hnw,h1,z1,h2,z2,b1,n1,c1,rho1,rhm1, &
!$omp                     x1,x2,x3,tauy,fd,ff,fb1,fb2,fb,taub)
do j=1,ni_vect
  i=i_vect(j)
  grad(i)=sign(min(abs(grad(i)),tanp),grad(i))
  cost=cos(atan(grad(i)))
  if(uo(i)>=0.d0)then ! when uo=0
    hnw=h(i)
    if(hnw<=0.d0)then
      u(i)=0.d0
      cycle
    end if
    h1=ho(i)
    z1=z(i)
    h2=ho(i_to(i))
    z2=z(i_to(i))
    b1=b(i)
    n1=n(i)
    c1=c(i)
    rho1=rho(i)
    rhm1=rhm(i)
    if(h1<hlim .or. z1+h1<z2)then
      u(i)=0.d0
      cycle
    end if
    where(i_from(i,:)==0)
      uhb(:)=0.d0
    elsewhere
      uhb(:)=uo(i_from(i,:))*ho(i_from(i,:))*b(i_from(i,:))
    end where
    x1=-uo(i)*(uo(i)*ho(i)*b(i)-sum(uhb(:)))/(dx(i)*b(i))
  else
    hnw=h(i_to(i))
    if(hnw<=0.d0)then
      u(i)=0.d0
      cycle
    end if
    h1=ho(i_to(i))
    h2=ho(i)
    z2=z(i)
    b1=b(i_to(i))
    n1=n(i_to(i))
    c1=c(i_to(i))
    rho1=rho(i_to(i))
    rhm1=rhm(i_to(i))
    if(h1<hlim .or. z1+h1<z2)then
      u(i)=0.d0
      cycle
    end if
    x1=-uo(i)*(uo(i_to(i))*h1*b1-uo(i)*ho(i)*b(i))/(dx(i)*b(i))
  end if
  !
  x2=-g*h1*grad(i)
  !
  if(c1<0.05d0)then
    !x3=-g*n1*n1*uo(i)*abs(uo(i))/h1**f1d3
    fb=(6.d0+2.5d0*log((h1/dm)))**(-2.d0)
  else
    fd=kd*sig/rho1*(1.d0-e2)*c1**(f1d3)
    ff=kf*(1.d0-c1)**f5d3/c1**(f2d3)
    fb1=6.25d0*(fd+ff)*((h1/dm))**(-2.d0)
    fb2=(6.d0+2.5d0*log((h1/dm)))**(-2.d0)
    fb=fb1
    if(fb1<fb2)fb=fb2
  end if
  tauy=(c1/csta)**f1d5*(sig-rho1)*c1*g*h1*cost*tanp
  taub=sign(tauy+rho1*fb*uo(i)*uo(i),uo(i))
  x3=-taub/rhm1
  u(i)=(uo(i)*h1-(-X1-X2+tauy/rhm1)*dt)/(hnw+rho1*fb*abs(uo(i))*dt/rhm1)
  !u(i)=(uo(i)*h1+(x1+x2+x3)*dt)/hnw
  if(u(i)*uo(i)<0.d0)then
    u(i)=0.d0
  end if
end do
!$omp end parallel do
umax=maxval(abs(u(:)))
iloc=maxloc(abs(u(:)))
i_umax=iloc(1)
humax=h(i_umax)
!
end subroutine cal_u
!
!
!
subroutine cal_q(ni_vect,i_vect,iend,i_to,b,h,u,q,cc,cf,qcc,qcf)
use omp_lib
implicit none
integer :: i,iend,j,ni_vect,i_vect(ni_vect),i_to(iend)
real(8) :: b(iend),h(iend),u(iend),q(0:iend),cc(iend),cf(iend),qcc(0:iend),qcf(0:iend)
!
!$omp parallel do private(i,j)
do j=1,ni_vect
  i=i_vect(j)
  if(u(i)>=0.d0)then
    q(i)=b(i)*h(i)*u(i)
    qcc(i)=q(i)*cc(i)
    qcf(i)=q(i)*cf(i)*(1.d0-cc(i))
  else
    q(i)=b(i_to(i))*h(i_to(i))*u(i)
    qcc(i)=q(i)*cc(i_to(i))
    qcf(i)=q(i)*cf(i_to(i))*(1.d0-cc(i_to(i)))
  end if
end do
!$omp end parallel do
!
end subroutine cal_q
!
!
!
subroutine cal_qs(ni_vect,i_vect,iend,i_to,b,hs,grad,qs,ks)
use omp_lib
implicit none
integer :: i,iend,j,ni_vect,i_vect(ni_vect),i_to(iend)
real(8) :: b(iend),hs(iend),grad(iend),qs(0:iend),ks,tant
!
!$omp parallel do private(i,j,tant)
do j=1,ni_vect
  i=i_vect(j)
  tant=grad(i)
  if(tant>0.d0)then
    qs(i)=sign(b(i_to(i))*ks*hs(i_to(i))*min(1.0d0,abs(tant)),-tant)
  else
    qs(i)=sign(b(i)*ks*hs(i)*min(1.0d0,abs(tant)),-tant)
  end if
end do
!$omp end parallel do
!
end subroutine cal_qs
!
!
!
subroutine cal_e(ni_scal,i_scal,iend,dt,dx,b,z,zs,sig,rho,csta,cstad,pf,pc,tanp,cc,cf,e,h,um,gradm, &
                 q_in,q_out,qcc_in,qcc_out,qcf_in,qcf_out,tant_e)
use omp_lib
implicit none
real(8), parameter :: PI=acos(-1.d0),deg2rad=PI/180.d0,rad2deg=180.d0/PI
integer :: i,iend,j,ni_scal,i_scal(ni_scal)
real(8) :: dt,sig,tanp,tant,tane,csta,cstad,pf,pc,emax,emin1,emin2,emin3, &
           z(iend),zs(iend),cc(iend),cf(iend),e(iend),h(iend),gradm(iend), &
           um(iend),rho(iend),tant_e(iend),q_in(iend),q_out(iend), &
           qcc_in(iend),qcc_out(iend),qcf_in(iend),qcf_out(iend),dx(iend),b(iend)
!
!$omp parallel do private(i,j,tant,tane,emax,emin1,emin2,emin3)
do j=1,ni_scal
  i=i_scal(j)
  if(um(i)>=0.d0)then
    tant=-gradm(i)
  else
    tant=gradm(i)
  end if
  tant=max(-gradm(i),0.d0)
  tant_e(i)=tant
  tane=(sig/rho(i)-1.d0)*cc(i)*tanp/((sig/rho(i)-1.d0)*cc(i)+1.d0)
  ! tan(a-b)=(tan(a)-tan(b))/(1+tan(a)*tan(b))
  e(i)=csta*(tant-tane)/(1.d0+tant*tane)*abs(um(i))
  if(e(i)>=0.d0)then ! ero
    emax=(z(i)-zs(i))/dt*csta*cos(atan(tant))
    e(i)=min(emax,e(i))
   ! if(ifld(i)==1)e(i)=emax
  else ! depo
    ! -h(i)/dt+(q(i+1)-q(i))/dx=e(i)/csta
    ! -cc(i)*h(i)/dt+(qcc(i+1)-qcc(i))/dx=e(i)
    ! -cf(i)*(1.d0-cc(i))*h(i)/dt+(qcf(i+1)-qcf(i))/dx=(1.d0-csta)/csta*cf(i)*e(i)
    emin1=min(0.d0,(-h(i)/dt*cstad+(q_out(i)-q_in(i))/dx(i)/b(i))*cstad)
    emin2=min(0.d0,-cc(i)*h(i)/dt+(qcc_out(i)-qcc_in(i))/dx(i)/b(i))
    if(cf(i)>0.d0)then
      emin3=min(0.d0,(-cf(i)*(1.d0-cc(i))*h(i)/dt+(qcf_out(i)-qcf_in(i))/dx(i)/b(i))*cstad/((1.d0-cstad)*cf(i)))
    else
      emin3=e(i)
    end if
    e(i)=max(emin1,emin2,emin3,e(i))
  end if
end do
!$omp end parallel do
!
end subroutine cal_e
!
!
!
subroutine cal_h(ni_scal,i_scal,iend,dt,dx,b,h,iqnum,q_in,q_out,e,csta,cstad)
use omp_lib
implicit none
integer :: i,iend,j,ni_scal,i_scal(ni_scal),iqnum
real(8) :: dt,csta,cstad,b(iend),h(iend),e(iend),dx(iend),q_in(iqnum),q_out(iqnum)
!
!$omp parallel do private(i,j)
do j=1,ni_scal
  i=i_scal(j)
  if(e(i)>=0.d0)then
    h(i)=h(i)-(q_out(i)-q_in(i))/dx(i)/b(i)*dt+e(i)/csta*dt
  else
    h(i)=h(i)-(q_out(i)-q_in(i))/dx(i)/b(i)*dt+e(i)/cstad*dt
  end if
  if(h(i)<0.d0)then
    if(h(i)<-0.001d0)then
      write(*,'(a,f10.3,i5,f12.7)')'h<0 @ (i)',h(i),i,e(i)
      write(*,*)q_out(i),q_in(i)
      !read(*,*)
    end if
    h(i)=0.d0
  end if
end do
!$omp end parallel do
!
end subroutine cal_h
!
!
!
subroutine cal_hs(ni_scal,i_scal,iend,dt,dx,rin,lam,b,hs,iqnum,qs_in,qs_out)
use omp_lib
implicit none
integer :: i,iend,j,ni_scal,i_scal(ni_scal),iqnum
real(8) :: dt,lam,b(iend),hs(iend),rin(iend),dx(iend),qs_in(iqnum),qs_out(iqnum)
!
!$omp parallel do private(i,j)
do j=1,ni_scal
  i=i_scal(j)
  hs(i)=hs(i)-(qs_out(i)-qs_in(i))/dx(i)/b(i)/lam*dt+rin(i)/lam*dt
  if(hs(i)<0.d0)then
    !write(*,'(a,i5,3f12.7)')'hs<0 @ (i)',i,hs(i),qs_out(i),qs_in(i)
    !read(*,*)
  end if
end do
!$omp end parallel do
!
end subroutine cal_hs
!
!
!
subroutine cal_sf(ni_scal,i_scal,iend,g,rhow,sig,tanp,lam,hlim,d,grad,hs,h,pw,c,sf,hsc)
use omp_lib
implicit none
integer :: i,iend,j,ni_scal,i_scal(ni_scal)
real(8) :: g,rhow,sig,lam,hlim,tant,cost,sint,tanp,ctmp,fg,fr,hsc0, &
           d(iend),grad(iend),hs(iend),h(iend),pw(iend),c(iend),sf(iend),hsc(iend)
!
!$omp parallel do private(i,j,tant,cost,sint,fg,fr,hsc0)
do j=1,ni_scal
  i=i_scal(j)
  tant=-grad(i)
  cost=cos(atan(tant))
  sint=sin(atan(tant))
  if(d(i)>=hlim)then
    ctmp=c(i)/(rhow*g*d(i)*cost*tanp)
    fg=rhow*g*d(i)*sint*(sig/rhow*(1.d0-lam)+(1.d0-hs(i)/d(i))*pw(i)+hs(i)/d(i)*lam+h(i)/d(i))
    fr=rhow*g*d(i)*cost*(sig/rhow*(1.d0-lam)+(1.d0-hs(i)/d(i))*pw(i)-hs(i)/d(i)*(1.d0-lam))*tanp+c(i)
    sf(i)=fg/fr
    !if(sf(i)>1.d0)ifld(i)=1
    hsc0=((1.d0-tant/tanp)*((1.d0-lam)*sig/rhow+pw(i))+ctmp) &
    /((1.d0-tant/tanp)*((1.d0-lam)+pw(i))+tant/tanp)
    hsc(i)=max(0.d0,min(1.d0,hsc0/d(i)))
  else
    sf(i)=1.d0
    hsc(i)=1.d0
  end if
end do
!$omp end parallel do
!
end subroutine cal_sf
!
!
!
subroutine cal_cc(ni_scal,i_scal,iend,dt,dx,b,h,ho,cc,cco,iqnum,q_in,q_out,e,gam,pc)
use omp_lib
implicit none
integer :: i,iend,j,ni_scal,i_scal(ni_scal),iqnum
real(8) :: dt,gam,pc,b(iend),h(iend),ho(iend),cc(iend),cco(iend),e(iend), &
           dx(iend),q_in(iqnum),q_out(iqnum)
!
!$omp parallel do private(i,j)
do j=1,ni_scal
  i=i_scal(j)
  if(h(i)<=0.d0)then ! checking h(i)==hmin
    cc(i)=0.d0
  else
    if(e(i)>=0.d0)then
      cc(i)=(cco(i)*ho(i)-(q_out(i)-q_in(i))*gam/dx(i)/b(i)*dt+pc*e(i)*dt)/h(i)
    else
      cc(i)=(cco(i)*ho(i)-(q_out(i)-q_in(i))*gam/dx(i)/b(i)*dt+e(i)*dt)/h(i)
    end if
  end if
  if(cc(i)<0.d0)then
    if(cc(i)*h(i)<-0.0001d0)then
      write(*,'(a,i5,5f15.7)')'cc*h<0 @ i,cc,h,e',i,cc(i),h(i),e(i)*dt/h(i)
      !read(*,*)
    end if
    cc(i)=0.d0
  end if
end do
!$omp end parallel do
!
end subroutine cal_cc
!
!
!
subroutine cal_cf(ni_scal,i_scal,iend,dt,dx,b,h,ho,cf,cfo,cc,cco,iqnum,q_in,q_out,e,gam,pf,cstad)
use omp_lib
implicit none
integer :: i,iend,j,ni_scal,i_scal(ni_scal),iqnum
real(8) :: dt,gam,pf,cstad,b(iend),h(iend),ho(iend),cf(iend),cfo(iend),e(iend), &
           cc(iend),cco(iend),dx(iend),q_in(iqnum),q_out(iqnum)
!
!$omp parallel do private(i,j)
do j=1,ni_scal
  i=i_scal(j)
  if(h(i)<=0.d0)then ! checking h(i)==hmin
    cf(i)=0.d0
  else
    if(e(i)>=0.d0)then
      cf(i)=(cfo(i)*(1.d0-cco(i))*ho(i) &
            -(q_out(i)-q_in(i))*gam/dx(i)/b(i)*dt+pf*e(i)*dt) &
             /((1.d0-cc(i))*h(i))
    else
      cf(i)=(cfo(i)*(1.d0-cco(i))*ho(i) &
            -(q_out(i)-q_in(i))*gam/dx(i)/b(i)*dt &
            +((1.d0-cstad)/cstad)*cfo(i)*e(i)*dt) &
            /((1.d0-cc(i))*h(i))
    end if
  end if

  if(cf(i)<0.d0)then
    if(cf(i)<-0.0001d0)then
      write(*,'(a,i5,5f15.7)')'cf<0 @ i,cf,h',i,cf(i),h(i)
      read(*,*)
    end if
    cf(i)=0.d0
  end if
end do
!$omp end parallel do
!
end subroutine cal_cf
!
!
!
subroutine cal_z(ni_scal,i_scal,iend,dt,csta,z,zs,e,tant_e)
use omp_lib
implicit none
real(8), parameter :: PI=acos(-1.d0),deg2rad=PI/180.d0,rad2deg=180.d0/PI
integer :: i,iend,j,ni_scal,i_scal(ni_scal)
real(8) :: dt,csta,tant,cost,z(iend),zs(iend),e(iend),tant_e(iend)
!
!$omp parallel do private(i,j,tant,cost)
do j=1,ni_scal
  i=i_scal(j)
  !tant=slope(i)
  tant=tant_e(i)
  cost=cos(atan(tant))
  z(i)=z(i)-e(i)/csta/cost*dt
  if(z(i)<zs(i))then
    !write(*,'(a,i5,5f10.3)')'z < zs at ',i,z(i),zs(i),e(i)
    !read(*,*)
    z(i)=zs(i)
  end if
end do
!$omp end parallel do
!
end subroutine cal_z
!
!
!
subroutine cal_vout(ni_scal,i_scal,iend,dt,dx,b,h,z1,z2,cc,cf,iacc,iaccmax,&
                    va_out,vc_out,vf_out)
use omp_lib
implicit none
real(8), parameter :: u2=0.35d0,g=9.8d0
integer :: ni_scal,i_scal(ni_scal),i,iend,j,iacc(iend),iaccmax
real(8) :: dt,wl,h1,qout,vtmp,b(iend),dx(iend),h(iend),z1(iend),z2(iend),&
           cc(iend),cf(iend),va_out(iend),vc_out(iend),vf_out(iend)
!
!$omp parallel do private(i,j,wl,h1,qout,vtmp)
do j=1,ni_scal
  i=i_scal(j)
  wl=z1(i)+h(i)
  if(z2(i)<wl .and. iaccmax<iacc(i))then
    h1=min(h(i),wl-z2(i))
    qout=u2*h1*(2.d0*g*h1)**0.5d0
    vtmp=qout*dt*2.d0*dx(i)
    va_out(i)=va_out(i)+vtmp
    vc_out(i)=vc_out(i)+vtmp*cc(i)
    vf_out(i)=vf_out(i)+vtmp*(1.d0-cc(i))*cf(i)
    h(i)=h(i)-vtmp/(b(i)*dx(i))
    if(h(i)<0.d0)then
      write(*,*)h(i),wl,z2(i),z1(i)
      read(*,*)
    end if
  end if
end do
!$omp end parallel do
!
end subroutine cal_vout
