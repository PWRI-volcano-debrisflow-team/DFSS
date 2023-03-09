! Program: RR_ver_1.0.f90
! Copyright (2021) by Public Works Research Institute (P.W.R.I.)
! License: CC-BY-SA
!-------------------------------------------------------------------------------
!     1    2    3      end  end+1
! FX  q    q    q       q   q
!     ->   ->  ->   ->  ->  ->
! CV  | z  | z  | z |...| z |
!     | 1  | 2  | 3 |...|end|
!-------------------------------------------------------------------------------
use omp_lib
implicit none
real(8), parameter :: PI=acos(-1.0d0),deg2rad=PI/180.0d0,rad2deg=180.0d0/PI, &
                      g=9.8d0,dtlim=0.0001d0,cfl=0.01d0
real(8), allocatable :: z(:,:),n(:,:),zs(:,:),d(:,:),lamda(:,:),pwc(:,:), &
                        kk(:,:),finf1(:,:),finf2(:,:),d_vash(:,:),facc(:,:)
real(8), allocatable :: h(:,:),ho(:,:),hs(:,:),hso(:,:),pw(:,:),pwo(:,:), &
                        eva1(:,:),eva2(:,:),ha(:,:),htmp(:,:),dqd(:,:),dqsd(:,:), &
                        hmax(:,:),hsmax(:,:),grdmax(:,:),grdave(:,:)
real(8), allocatable :: qx(:,:),qy(:,:),qxo(:,:),qyo(:,:), &
                        qsx(:,:),qsy(:,:),qsxo(:,:),qsyo(:,:), &
                        qx_acc(:,:),qy_acc(:,:),qsx_acc(:,:),qsy_acc(:,:), &
                        umax(:,:),vmax(:,:),h_umax(:,:),h_vmax(:,:)
integer, allocatable :: i_2d(:,:),j_2d(:,:),iriver(:,:),ibasin(:,:), &
                        iacc(:,:),idrain(:,:)
logical, allocatable :: m_cv(:,:)
integer :: n_ij_cv,n_ij_u,n_ij_v,n_ij_u_we,n_ij_v_sn
integer, allocatable :: ij_cv(:,:),ij_u(:,:),ij_v(:,:),ij_u_we(:,:),ij_v_sn(:,:)
real(8), allocatable :: rseries(:,:,:),r(:,:),rmap(:,:,:),x_rpnt(:),y_rpnt(:), &
                        r_rpnt(:),d_idw(:),eva_tmp(:)
integer :: n_rpnt,minl(1)
character(len=100) :: rname
real(8) :: dx,dy,dt,rdt,output_f_time,output_d_time,xllcorner,yllcorner,cellsize, &
           vad0,d_uni,rn_uni,fs1,fs2,lam_uni,k_uni,pwc_uni,hs0,pw0,hthr
real(8) :: time,time0,time1,time2,end_time,dt0,dtdiv,dtmin,dtmin_h2d,flxd, &
           umax_disp,vmax_disp,h_umax_disp,h_vmax_disp, &
           r_total,qs_total,q_total,eva1_total,eva2_total,out_total,out_total_ini,finf2_total, &
           x,y,dist,deg,basin_area,rain_tmp,tmp1
integer :: iend,jend,itsp_max,itsp_f_out,itsp_d_out,itsp_rain,rnum
integer :: i,j,k,itsp,itr,itr2,maxloc_2d(2),i_umax,j_umax,i_vmax,j_vmax
character(len=100) :: fname,cfmt,ctmp,fotime,out_dir
character(len=18),allocatable :: datetime(:)
! ---> 1d
real(8), allocatable :: zini_1d(:),z_1d(:),zs_1d(:),b_1d(:),d_1d(:),n_1d(:), &
                        h_1d(:),ho_1d(:),q_in(:),q_out(:),dx_scal(:)
real(8), allocatable :: u_1d(:),uo_1d(:),q_1d(:),qo_1d(:),grad_z_1d(:),grad_w_1d(:),dx_vect(:)
integer :: ni_inlet,ni_outlet
integer, allocatable :: i_inlet(:),i_outlet(:)
!real(8), allocatable :: qin0(:)
integer :: ni_scal,ni_vect
integer, allocatable :: i_node(:),i_to(:),i_from(:,:),i_scal(:),i_vect(:)
integer, allocatable :: ij_2d1d(:,:),iacc_1d(:)
logical, allocatable :: mask_1d(:)
real(8) :: rn_river,wid_m,wid_p,wid_min,dep_m,dep_p,dep_min,hlim_1d
real(8) :: umax_1d,dtmin_h1d,h1,h2,ba,q_1d_total
real(8), allocatable :: x_2d_1d(:),y_2d_1d(:),vsl_out(:)
integer :: iend_1d,iflg_1d
! <-- CHECK Used CPU No.
! ref https://tama.green.gifu-u.ac.jp/~tama/Memo/openmp.html
 write(*,*) 'CHECK HOW MANY CPU USED'
 write(*,*) 'before parallel'
 !$omp parallel
 write(*,*) 'Hello OpenMP'
 !$omp end parallel
 write(*,*) 'after parallel'
 write(*,*) ''
 write(*,*) 'START'
! ---> CPU
! <--- 1d
!
open(1,file='RR_input.txt')
read(1,*)
read(1,*)iend
read(1,*)jend
read(1,*)dx
read(1,*)dy
!
allocate(z(iend,jend),n(iend,jend),zs(iend,jend),d(iend,jend), &
         lamda(iend,jend),pwc(iend,jend),kk(iend,jend),finf1(iend,jend), &
         finf2(iend,jend),d_vash(iend,jend),facc(iend,jend))
allocate(h(iend,jend),ho(iend,jend),hs(iend,jend),hso(iend,jend), &
         pw(iend,jend),pwo(iend,jend),eva1(iend,jend),eva2(iend,jend), &
         ha(iend,jend),htmp(iend,jend),dqd(iend,jend),dqsd(iend,jend), &
         hmax(iend,jend),hsmax(iend,jend),grdmax(iend,jend),grdave(iend,jend))
allocate(qx(iend+1,jend+1),qy(iend+1,jend+1),qxo(iend+1,jend+1),qyo(iend+1,jend+1), &
         qsx(iend+1,jend+1),qsy(iend+1,jend+1),qsxo(iend+1,jend+1),qsyo(iend+1,jend+1), &
         qx_acc(iend+1,jend+1),qy_acc(iend+1,jend+1),qsx_acc(iend+1,jend+1),qsy_acc(iend+1,jend+1), &
         umax(iend+1,jend+1),vmax(iend+1,jend+1),h_umax(iend+1,jend+1),h_vmax(iend+1,jend+1))
allocate(i_2d(iend+1,jend+1),j_2d(iend+1,jend+1),iriver(iend,jend), &
         ibasin(iend,jend),iacc(iend,jend),idrain(iend,jend),m_cv(iend,jend))
!
read(1,*)
read(1,*)dt
read(1,*)output_f_time
read(1,*)output_d_time
read(1,*)
read(1,*)rname
read(1,*)rnum
read(1,*)rdt
!
allocate(rseries(iend,jend,0:rnum))
allocate(eva_tmp(0:rnum),datetime(0:rnum))
!
read(1,*)
read(1,*)fname
write(*,*)'elevation',trim(adjustl(fname))
call r_fasc(10,fname,z,iend,jend,xllcorner,yllcorner,cellsize)
!
read(1,*)fname
write(*,*)'targetArea',trim(adjustl(fname))
call r_iasc(10,fname,ibasin,iend,jend,xllcorner,yllcorner,cellsize)
basin_area=float(count(ibasin(:,:)/=0))*dx*dx
write(*,'(a,f15.3)')' area (m2)',basin_area
!
read(1,*)fname
write(*,*)'flowAcc',trim(adjustl(fname))
call r_fasc(10,fname,facc,iend,jend,xllcorner,yllcorner,cellsize)

read(1,*)fname
write(*,*)'flowDir',trim(adjustl(fname))
call r_iasc(10,fname,idrain,iend,jend,xllcorner,yllcorner,cellsize)
!
read(1,*)fname
write(*,*)'volcanicAsh',trim(adjustl(fname))
call r_fasc(10,fname,d_vash,iend,jend,xllcorner,yllcorner,cellsize)
read(1,*)vad0 !depth of ash(m) when fs=0
!
read(1,*)
read(1,*)i,d_uni,fname !soil depth(m)
if(i/=0)then
  write(*,*)'d         <--- ',trim(adjustl(fname))
  call r_fasc(10,fname,d,iend,jend,xllcorner,yllcorner,cellsize)
else
  d(:,:)=d_uni
end if
fname=trim(adjustl('output_RR_misc/depth.asc'))
cfmt='(*(f10.3))'
call w_fasc(20,fname,cfmt,d,iend,jend,xllcorner,yllcorner,dx,-9999.0d0)
!
read(1,*)i,rn_uni,fname !Maning's roughness
if(i/=0)then
  write(*,*)'n         <--- ',trim(adjustl(fname))
  call r_fasc(10,fname,n,iend,jend,xllcorner,yllcorner,cellsize)
else
  n(:,:)=rn_uni
end if
fname=trim(adjustl('output_RR_misc/roughness.asc'))
cfmt='(*(f10.3))'
call w_fasc(20,fname,cfmt,n,iend,jend,xllcorner,yllcorner,dx,-9999.0d0)
!
read(1,*)i,fs1,fname !surface infiltration(mm/h)
if(i/=0)then
  write(*,*)'fs1(cm/s) <--- ',trim(adjustl(fname))
  call r_fasc(10,fname,finf1,iend,jend,xllcorner,yllcorner,cellsize)
  finf1(:,:)=finf1(:,:)/3600.0d0/1000.0d0 !change unit (mm/h) -> (m/s)
else
  finf1(:,:)=fs1/3600.0d0/1000.0d0        !change unit (mm/h) -> (m/s)
end if
!
read(1,*)i,fs2,fname ! infiltration to lower layer(mm/h)
if(i/=0)then
  write(*,*)'fs2(cm/s) <--- ',trim(adjustl(fname))
  call r_fasc(10,fname,finf2,iend,jend,xllcorner,yllcorner,cellsize)
  finf2(:,:)=finf2(:,:)/3600.0d0/1000.0d0 !change unit (mm/h) -> (m/s)
else
  finf2(:,:)=fs2/3600.0d0/1000.0d0        !change unit (mm/h) -> (m/s)
end if
!
finf1(:,:)=(-d_vash(:,:)/vad0+1.d0)*finf1(:,:)
where(finf1(:,:)<0.d0)finf1(:,:)=0.d0
where(ibasin(:,:)/=0)finf2(:,:)=0.d0 !fs2
fname=trim(adjustl('output_RR_misc/fs1_mm_h.asc'))
cfmt='(*(f10.3))'
call w_fasc(20,fname,cfmt,finf1*3.6d+6,iend,jend,xllcorner,yllcorner,dx,-9999.0d0)
fname=trim(adjustl('output_RR_misc/fs2_mm_h.asc'))
cfmt='(*(f10.3))'
call w_fasc(20,fname,cfmt,finf2*3.6d+6,iend,jend,xllcorner,yllcorner,dx,-9999.0d0)
!
read(1,*)i,lam_uni,fname !porosity
if(i/=0)then
  write(*,*)'lamda     <--- ',trim(adjustl(fname))
  call r_fasc(10,fname,lamda,iend,jend,xllcorner,yllcorner,cellsize)
else
  lamda(:,:)=lam_uni
end if
fname=trim(adjustl('output_RR_misc/lamda.asc'))
cfmt='(*(f10.3))'
call w_fasc(20,fname,cfmt,lamda,iend,jend,xllcorner,yllcorner,dx,-9999.0d0)
!
read(1,*)i,k_uni,fname !saturated hydraulic conductivity(cm/s)
if(i/=0)then
  write(*,*)'k(cm/s)   <--- ',trim(adjustl(fname))
  call r_fasc(10,fname,kk,iend,jend,xllcorner,yllcorner,cellsize)
else
  kk(:,:)=k_uni
end if
kk=kk/100.d0
fname=trim(adjustl('output_RR_misc/k_cm_s.asc'))
cfmt='(*(e10.3))'
call w_fasc(20,fname,cfmt,kk*1.d+2,iend,jend,xllcorner,yllcorner,dx,-9999.0d0)
!
read(1,*)i,pwc_uni,fname !clitic vwc for subsurface flow initiation
if(i/=0)then
  write(*,*)'pwc       <--- ',trim(adjustl(fname))
  call r_fasc(10,fname,pwc,iend,jend,xllcorner,yllcorner,cellsize)
else
  pwc(:,:)=pwc_uni
end if
fname=trim(adjustl('output_RR_misc/pwc.asc'))
cfmt='(*(f10.3))'
call w_fasc(20,fname,cfmt,pwc,iend,jend,xllcorner,yllcorner,dx,-9999.0d0)
!
read(1,*)i,pw0,fname !initial soil vwc
if(i/=0)then
  write(*,*)'pw_ini    <--- ',trim(adjustl(fname))
  call r_fasc(10,fname,pw,iend,jend,xllcorner,yllcorner,cellsize)
else
  pw(:,:)=pw0
end if
where(ibasin(:,:)==0)pw(:,:)=0.d0
fname=trim(adjustl('output_RR_misc/pw_ini.asc'))
cfmt='(*(f10.3))'
call w_fasc(20,fname,cfmt,pw,iend,jend,xllcorner,yllcorner,dx,-9999.0d0)
!
read(1,*)i,hs0,fname !initial hs(m)
if(i/=0)then
  write(*,*)'hs_ini    <--- ',trim(adjustl(fname))
  call r_fasc(10,fname,hs,iend,jend,xllcorner,yllcorner,cellsize)
  else
    hs=hs0
end if
where(ibasin(:,:)==0)hs(:,:)=0.d0
fname=trim(adjustl('output_RR_misc/hs_ini.asc'))
cfmt='(*(f10.3))'
call w_fasc(20,fname,cfmt,hs,iend,jend,xllcorner,yllcorner,dx,-9999.0d0)
!
zs(:,:)=z(:,:)-d(:,:)
call cal_grad(iend,jend,dx,zs,grdave)
fname=trim(adjustl('output_RR_misc/gradave.asc'))
cfmt='(*(f10.3))'
call w_fasc(20,fname,cfmt,grdave*rad2deg,iend,jend,xllcorner,yllcorner,dx,-9999.0d0)
!
read(1,*)
read(1,*)hthr !threshold of water depth to move(m)
!
! ---> stream
read(1,*)
read(1,*)iflg_1d
read(1,*)rn_river
read(1,*)wid_m
read(1,*)wid_p
read(1,*)wid_min
read(1,*)dep_m
read(1,*)dep_p
read(1,*)dep_min
! <--- stream
close(1)
write(*,*)'ctl read end'
!
end_time=dble(rnum-1)*rdt
itsp_max=nint(end_time/dt)
itsp_f_out=nint(output_f_time/dt)
itsp_d_out=nint(output_d_time/dt)
itsp_rain=nint(rdt/dt)
!
out_dir='output_RR/'
write(*,*)trim(adjustl(out_dir))
write(*,*)''
! ---> control points etc.
do i=1,iend+1
  do j=1,jend+1
    i_2d(i,j)=i
    j_2d(i,j)=j
  end do
end do
!
m_cv(:,:)=(ibasin(:,:)/=0)
!
open(1,file='controlVolume_CenterPoint_inRR.txt')
read(1,*)n_ij_cv
allocate(ij_cv(n_ij_cv,2))
read(1,*)
do k=1,n_ij_cv
  read(1,*)x,y,ij_cv(k,1:2)
end do
close(1)
open(1,file='fluxPoint_x_inRR.txt')
read(1,*)n_ij_u
allocate(ij_u(n_ij_u,2))
read(1,*)
do k=1,n_ij_u
  read(1,*)x,y,ij_u(k,1:2)
end do
close(1)
open(1,file='fluxPoint_y_inRR.txt')
read(1,*)n_ij_v
allocate(ij_v(n_ij_v,2))
read(1,*)
do k=1,n_ij_v
  read(1,*)x,y,ij_v(k,1:2)
end do
close(1)
open(1,file='fluxPoint_xBoudary_inRR.txt')
read(1,*)n_ij_u_we
allocate(ij_u_we(n_ij_u_we,8))
read(1,*)
do k=1,n_ij_u_we
  read(1,*)x,y,ij_u_we(k,:)
end do
close(1)
open(1,file='fluxPoint_yBoudary_inRR.txt')
read(1,*)n_ij_v_sn
allocate(ij_v_sn(n_ij_v_sn,8))
read(1,*)
do k=1,n_ij_v_sn
  read(1,*)x,y,ij_v_sn(k,:)
end do
close(1)
! <--- control points etc.
!
! ---> rain and eva
open(11,file=rname)
read(11,*)n_rpnt
allocate(r(iend,jend),rmap(iend,jend,0:n_rpnt),x_rpnt(n_rpnt),y_rpnt(n_rpnt), &
         r_rpnt(n_rpnt),d_idw(n_rpnt))
read(11,*)ctmp,x_rpnt(1:n_rpnt)
read(11,*)ctmp,y_rpnt(1:n_rpnt)
rseries(:,:,0:rnum)=0.0d0 ! from 0
rmap(:,:,:)=0
do i=1,iend
  do j=1,jend
    if(ibasin(i,j)/=0)then
      x=xllcorner+dx*(0.5+i-1)
      y=yllcorner+dx*(0.5+j-1)
      d_idw(1:n_rpnt)=((x-x_rpnt(1:n_rpnt))**2.d0+(y-y_rpnt(1:n_rpnt))**2.d0)**0.5d0
      minl=minloc(d_idw)
      !write(*,*)i,j,minl
      rmap(i,j,1:n_rpnt)=0.d0
      rmap(i,j,minl(1))=1.d0 !d_idw(1:n_rpnt)/sum(d_idw(:))
      rmap(i,j,0)=minl(1)
    end if
  end do
end do
fname='output_RR_misc/rain_map.asc'
cfmt='(*(f10.3))'
call w_fasc(20,fname,cfmt,rmap(:,:,0),iend,jend,xllcorner,yllcorner,dx,-9999.d0)
!
do k=1,n_rpnt
  write(ctmp,'(i4.4)')k
  fname='output_RR_misc/rain_map'//trim(ctmp)//'.asc'
  cfmt='(*(f10.3))'
  !call w_fasc(20,fname,cfmt,rmap(:,:,k),iend,jend,xllcorner,yllcorner,dx,-9999.d0)
end do
!
open(20,file='output_RR_misc/basin_rain.txt')
write(20,'(a)')'datetime,rain(mm),maxval'
!
rseries(:,:,:)=0.0d0
do k=0,rnum-1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  read(11,*)datetime(k),r_rpnt(1:n_rpnt),eva_tmp(k)
  if(k==0 .or. k==rnum)write(*,'(x,a,a,i5)')datetime(k),' rain # ',k
  eva_tmp(k)=max(eva_tmp(k),0.d0)/rdt/1000.0d0
    do i=1,iend
      do j=1,jend
        if(ibasin(i,j)/=0)then
          !rain_tmp=sum(r_rpnt(1:n_rpnt)*rmap(i,j,1:n_rpnt))
          rain_tmp=r_rpnt(int(rmap(i,j,0)))
          !write(*,'(2i4,10f7.3)')i,j,rain_tmp,r_rpnt(:),rmap(i,j,:)
          rseries(i,j,k)=rain_tmp/rdt/1000.0d0
        end if
      end do
    end do
    !write(*,*)datetime(k),maxval(rseries(:,:,k),mask=ibasin/=0)*rdt*1000.d0
    write(20,'(a,2(a,f10.3))')trim(adjustl(datetime(k))),',', &
    sum(rseries(:,:,k),mask=ibasin/=0)*dx*dx/basin_area*rdt*1000.0d0,',',maxval(rseries(:,:,k),mask=ibasin/=0)*rdt*1000.d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!rseries(:,:,k)=rain_tmp/rdt/1000.0d0
!where(ibasin(:,:)==0) rseries(:,:,k)=0.0d0
!write(fname,'(i10.10)')nint(k*rdt)
!fname='output_rain/rain_'//trim(adjustl(fname))//'.asc'
!write(*,*)trim(adjustl(fname)),' mm/h'
!cfmt='(*(f10.3))'
!call w_fasc(20,cfmt,rseries(:,:,k)*3600.*1000.0d0,iend,jend,xllcorner,yllcorner,dx,-9999.0d0)
!
!open(10,file=fname)
!do i=1,6
!  read(10,*)
!end do
!do j=jend,1,-1
!  read(10,*)rseries(1:iend,j,k)
!  rseries(1:iend,j,k)=rseries(1:iend,j,k)/rdt/1000.0d0
!end do
!
end do
close(11)
close(20)
!
fname='output_RR_misc/rain_sum.asc'
r=0.d0
do k=0,rnum-1
 r(:,:)=r(:,:)+(rseries(:,:,k))*rdt*1000.
end do
write(*,'(a,f10.3,a)')' basin ave. rain: ',sum(r(:,:))*dx*dx/basin_area,' mm'
write(*,*)''
cfmt='(*(f10.3))'
call w_fasc(20,fname,cfmt,r(:,:),iend,jend,xllcorner,yllcorner,dx,-9999.0d0)
r=0.d0
rseries(:,:,rnum)=rseries(:,:,rnum-1)
datetime(rnum)=datetime(rnum-1)
eva_tmp(rnum)=eva_tmp(rnum-1)
! <--- rain and eva
!
! ---> 1d
open(10,file='streamConfiguration_inRR.txt')
read(10,*)iend_1d
i=iend_1d
!
allocate(zini_1d(i),z_1d(i),zs_1d(i),b_1d(i),d_1d(i),n_1d(i),h_1d(i),ho_1d(i), &
         q_in(i),q_out(i),iacc_1d(i),dx_scal(i))
allocate(u_1d(i),uo_1d(i),q_1d(0:i),qo_1d(0:i),grad_z_1d(i),grad_w_1d(i),dx_vect(i))
!allocate(qin0(0:qnum),cin0(0:qnum))
allocate(i_node(i),i_to(i),i_from(i,3),mask_1d(i))
allocate(ij_2d1d(i,2),x_2d_1d(i),y_2d_1d(i),vsl_out(i))
!
read(10,*)
do i=1,iend_1d
  read(10,*)i_node(i),i_to(i),i_from(i,1),i_from(i,2),i_from(i,3), &
            tmp1,z_1d(i),zs_1d(i),dx_scal(i),b_1d(i),n_1d(i),iacc_1d(i), &
            x_2d_1d(i),y_2d_1d(i),ij_2d1d(i,1),ij_2d1d(i,2)
  dx_vect(i)=dx_scal(i)
  tmp1=z(ij_2d1d(i,1),ij_2d1d(i,2))
  if(iflg_1d==0)then
    ba=cellsize*cellsize*float(iacc_1d(i))/1000000.d0
    b_1d(i)=max(wid_min,wid_m*ba**wid_p)
    d_1d(i)=max(dep_min,dep_m*ba**dep_p)
    z_1d(i)=tmp1-d_1d(i)
    n_1d(i)=rn_river
  end if
end do
close(10)
!
zini_1d(:)=z_1d(:)
mask_1d=(i_from(:,1)/=0)! .and. i_to(:)/=0)
ni_scal=count(mask_1d)
allocate(i_scal(ni_scal))
i_scal(1:ni_scal)=pack(i_node(:),mask_1d)
open(20,file='output_RR_misc/scal_point.txt')
write(20,*)'x(s(i)),y(s(i)),s(i),i'
do i=1,ni_scal
  write(20,'(2(f15.5,:'',''),2(i5,:'',''))')x_2d_1d(i_scal(i)),y_2d_1d(i_scal(i)),i_scal(i),i
end do
close(20)
!
mask_1d=i_to(:)/=0
ni_vect=count(mask_1d)
allocate(i_vect(ni_vect))
i_vect(1:ni_vect)=pack(i_node(:),mask_1d)
open(20,file='output_RR_misc/vect_point.txt')
write(20,*)'x(v(i)),y(v(i)),dist,deg,i,i_to'
do i=1,ni_vect
  x=x_2d_1d(i_vect(i))
  y=y_2d_1d(i_vect(i))
  x=(x+x_2d_1d(i_to(i_vect(i))))*0.5d0
  y=(y+y_2d_1d(i_to(i_vect(i))))*0.5d0
  h1=x_2d_1d(i_to(i_vect(i)))-x_2d_1d(i_vect(i)) !dx->h1
  h2=y_2d_1d(i_to(i_vect(i)))-y_2d_1d(i_vect(i)) !dy->h2
  dist=(h1*h1+h2*h2)**0.5d0
  h1=h1/dist !dx->h1
  h2=h2/dist !dy->h2
  if(0<h1)then
    deg=atan(h2/h1)*rad2deg
  else if(h1<0)then
    deg=atan(h2/h1)*rad2deg
  else
    deg=90.
  end if
  deg=deg*(-45.d0)+90.d0
  write(20,'(2(f15.5,:'','')2(f10.3,:'',''),2(i5,:'',''))')&
        x,y,dist,deg,i_vect(i),i_to(i_vect(i))
end do
close(20)
!
mask_1d=i_from(:,1)==0
ni_inlet=count(mask_1d)
allocate(i_inlet(ni_inlet))
i_inlet(1:ni_inlet)=pack(i_node(:),mask_1d)
open(20,file='output_RR_misc/inlet_point.txt')
write(20,*)'x(s(i)),y(s(i)),s(i),i'
do i=1,ni_inlet
  write(20,'(2(f15.5,:'',''),2(i5,:'',''))')x_2d_1d(i_inlet(i)),y_2d_1d(i_inlet(i)),i_inlet(i),i
end do
close(20)
!
mask_1d=i_to(:)==0
ni_outlet=count(mask_1d)
allocate(i_outlet(ni_outlet))
i_outlet(1:ni_outlet)=pack(i_node(:),mask_1d)
open(20,file='output_RR_misc/outlet_point.txt')
write(20,*)'x(s(i)),y(s(i)),s(i),i'
do i=1,ni_outlet
  write(20,'(2(f15.5,:'',''),2(i5,:'',''))')x_2d_1d(i_outlet(i)),y_2d_1d(i_outlet(i)),i_outlet(i),i
end do
close(20)
!
u_1d(:)=0.d0
uo_1d(:)=u_1d(:)
h_1d(:)=0.d0
ho_1d(:)=h_1d(:)
q_1d(:)=0.d0
qo_1d(:)=q_1d(:)
q_in(:)=0.d0
q_out(:)=0.d0
hlim_1d=hthr
! <--- 1d
!
! ---> initial conditions
h(:,:)=0.d0
ho(:,:)=h(:,:)
hso(:,:)=hs(:,:)
hmax(:,:)=0.d0
hsmax(:,:)=0.d0
pwo(:,:)=pw(:,:)
eva1(:,:)=0.d0
eva2(:,:)=0.d0
!
qx(:,:)=0.0d0
qy(:,:)=0.0d0
qxo(:,:)=qx(:,:)
qyo(:,:)=qy(:,:)
qsx(:,:)=0.0d0
qsy(:,:)=0.0d0
qsxo(:,:)=qsx(:,:)
qsyo(:,:)=qsy(:,:)
qx_acc(1:iend+1,1:jend+1)=0.0d0
qy_acc(1:iend+1,1:jend+1)=0.0d0
qsx_acc(1:iend+1,1:jend+1)=0.0d0
qsy_acc(1:iend+1,1:jend+1)=0.0d0
!
r_total=0.0d0
q_total=0.0d0
qs_total=0.d0
eva1_total=0.d0
eva2_total=0.d0
finf2_total=0.d0
q_1d_total=0.d0
out_total_ini=finf2_total+qs_total+q_total+eva1_total+eva2_total+(sum(h)+sum(hs*(lamda-pwc))+sum(pw*d))*dx*dx
out_total=0.d0
vsl_out(:)=0.d0
umax_disp=0.d0
vmax_disp=0.d0
h_umax_disp=0.d0
h_vmax_disp=0.d0
maxloc_2d=0
i_umax=0
j_umax=0
i_vmax=0
j_vmax=0
dtmin=dt
dqd=0.d0
dqsd=0.d0
! <--- initial conditions
!
dt0=dt
dtmin=dt
time=0.0d0
time0=0.d0
open(200,file='RR_summary.txt')
write(200,'(10a21,4a11,a10)')&
'time,', 'rain,', 'out_t,', 'out_qs,', 'out_q,', 'out_eva1,', 'out_eva2,', &
'sum_hs,', 'sum_h,', 'sum_pw,', 'slope,', 'debris,','fs2,','out_q1d','out_h1d'
open(201,file='flowVolume_RRtoDR.txt')
write(201,'(a10,"",a10,*(i12,:,""))')'','i_2d',ij_2d1d(:,1)
write(201,'(a10,"",a10,*(i12,:,""))')'','j_2d',ij_2d1d(:,2)
write(200,'(a,i10,9(a,E20.10e3),2(a,i10),6(a,E20.10e3))')&
'  time(s)=',nint(time),',',r_total,',',out_total, &
',',qs_total,',',q_total,',',eva1_total,',',eva2_total,',', &
sum(hs*(lamda-pwc))*dx*dx,',',sum(h)*dx*dx,',',sum(pw*d)*dx*dx,',', &
0,',',0, ',',finf2_total,',',0.d0,',',0.d0
write(201,'(a10,"",i10,*(e12.3,:,""))')'time(s)=',nint(time),vsl_out(:)
!
write(*,*)'t=0 ---> ',datetime(0)
write(*,*)''
!
do itsp=1,itsp_max
  time=time0+dt
  k=itsp/itsp_rain
  ! --->  rain
  !r(:,:)=rseries(:,:,k)+(rseries(:,:,k+1)-rseries(:,:,k))*(itsp-k*itsp_rain)/dble(itsp_rain)
  r(:,:)=rseries(:,:,k) ! r=0,1,2,,,
  where(ibasin(:,:)==0) r(:,:)=0.d0
  ! <--- rain
  ! ---> eva
  eva1(:,:)=eva_tmp(k)+(eva_tmp(k+1)-eva_tmp(k))*(itsp-k*itsp_rain)/dble(itsp_rain)
  where(ibasin(:,:)==0) eva1(:,:)=0.d0
  eva2(:,:)=eva1(:,:)
  ! <--- eva
  time1=time0
  itr=1
  dt=dt0
  umax(:,:)=0.d0
  vmax(:,:)=0.d0
  h_umax(:,:)=0.d0
  h_vmax(:,:)=0.d0
  flxd=0.d0
  call cal_q(n_ij_u,ij_u,iend,jend,z,h+r*dt,n,-1,0,hthr,dx,qx,umax,h_umax,dt,0.d0)
  call cal_q(n_ij_v,ij_v,iend,jend,z,h+r*dt,n,0,-1,hthr,dx,qy,vmax,h_vmax,dt,0.d0)
  flxd=max(flxd,maxval(umax)*dt,maxval(vmax)*dt)
  call cal_qs(n_ij_u,ij_u,iend,jend,zs,hs,-1,0,kk,hthr,dx,qsx,umax,h_umax)
  call cal_qs(n_ij_v,ij_v,iend,jend,zs,hs,0,-1,kk,hthr,dx,qsy,vmax,h_vmax)
  flxd=max(flxd,maxval(umax)*dt,maxval(vmax)*dt)
  dtmin_h2d=dt
  call cal_dq(n_ij_cv,ij_cv,iend,jend,dx,h,qx,qy,dqd,dtmin_h2d) ! h as dq
  if(dt < dtmin_h2d)then
    write(*,*)'dt < dtmin_h2d',dt,dtmin_h2d
    dt=dtmin_h2d
    read(*,*)
  end if
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if ((dx*cfl) < flxd .or. dtmin_h2d < dt)then
    itr=max(1,ceiling(flxd/(dx*cfl)),nint(dt0/dtmin_h2d))
    if(itr<=0)then
      write(*,*)'itr',itr,flxd,dx*cfl,nint(dt0/dtlim)
      read(*,*)
    end if
    dt=dt0/dble(itr)
    write(*,'(a,f15.7,a,i10)', advance='yes')'2d dt=',dt,', itr=',itr
    if(dt<dtmin)dtmin=dt
  end if
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do while (itr>0)
    ! ---> cal 2D
    call cal_dq(n_ij_cv,ij_cv,iend,jend,dx,h,  qx, qy,dqd,0.d0)  ! h as dq
    call cal_dq(n_ij_cv,ij_cv,iend,jend,dx,hs,qsx,qsy,dqsd,0.d0) ! hs as dqs
    where(m_cv(:,:)) htmp(:,:)=min(ho(:,:)+(r(:,:)+dqd(:,:))*dt,finf1(:,:)*dt)
    where(m_cv(:,:)) ha(:,:)=htmp(:,:)+(lamda(:,:)-pwc(:,:))*hs(:,:)+d(:,:)*pw(:,:)+&
    !dt*(dqsd(:,:)-eva1(:,:)-eva2(:,:)) ! dt!
    dt*(dqsd(:,:)) ! dt!
    !where(m_cv(:,:))finf2(:,:)=fs2
    where(m_cv(:,:))finf2(:,:)=fs2*hs(:,:)
    call cal_divha(n_ij_cv,ij_cv,iend,jend,lamda,pwc,ha,h,hs,pw,d,finf2,dt,qx,qy,qsx,qsy)
    where(m_cv(:,:))htmp(:,:)=max(ho(:,:)+(r(:,:)+dqd(:,:))*dt-finf1(:,:)*dt,0.d0)
    h(:,:)=h(:,:)+htmp(:,:)
    ! ---> cal eva
    call cal_eva(n_ij_cv,ij_cv,iend,jend,r,hs,pw,lamda,pwc,d,dt,eva1,eva2)
    ! <--- cal eva
    !
    ! ---> cal flow from 2d to 1d
    call cal_2d1d(iend,jend,ij_2d1d,iend_1d,dx,dt,z,h,z_1d,h_1d,b_1d,dx_scal,vsl_out)
    ! <--- cal flow from 2d to 1d
    !
    ! ---> cal surface flux
    call cal_q(n_ij_u,ij_u,iend,jend,z,h,n,-1,0,hthr,dx,qx,umax,h_umax,dt,cfl)
    call cal_q(n_ij_v,ij_v,iend,jend,z,h,n,0,-1,hthr,dx,qy,vmax,h_vmax,dt,cfl)
    ! <--- cal surface flux
    !
    ! ---> cal subsurface flux
    call cal_qs(n_ij_u,ij_u,iend,jend,zs,hs,-1,0,kk,hthr,dx,qsx,umax,h_umax)
    call cal_qs(n_ij_v,ij_v,iend,jend,zs,hs,0,-1,kk,hthr,dx,qsy,vmax,h_vmax)
    ! <--- cal subsurface flux
    !
    qx_acc( 1:iend+1,1:jend+1)=qx_acc( 1:iend+1,1:jend+1)+qx( 1:iend+1,1:jend+1)*dx*dt
    qy_acc( 1:iend+1,1:jend+1)=qy_acc( 1:iend+1,1:jend+1)+qy( 1:iend+1,1:jend+1)*dx*dt
    qsx_acc(1:iend+1,1:jend+1)=qsx_acc(1:iend+1,1:jend+1)+qsx(1:iend+1,1:jend+1)*dx*dt
    qsy_acc(1:iend+1,1:jend+1)=qsy_acc(1:iend+1,1:jend+1)+qsy(1:iend+1,1:jend+1)*dx*dt
    ho(:,:)=h(:,:)
    hso(:,:)=hs(:,:)
    pwo(:,:)=pw(:,:)
    qxo(:,:)=qx(:,:)
    qyo(:,:)=qy(:,:)
    qsxo(:,:)=qsx(:,:)
    qsyo(:,:)=qsy(:,:)
    where(hmax(:,:)<h(:,:))hmax(:,:)=h(:,:)
    where(hsmax(:,:)<hs(:,:))hsmax(:,:)=hs(:,:)
    eva1_total=eva1_total+sum(eva1)*dx*dx*dt
    eva2_total=eva2_total+sum(eva2)*dx*dx*dt
    finf2_total=finf2_total+sum((pack(finf2(:,:),mask=m_cv(:,:)))*dt)*dx*dx
    do k=1,n_ij_u_we
      qs_total=qs_total+abs(-qsx(ij_u_we(k,1),ij_u_we(k,2)))*dx*dt
      q_total =q_total +abs(- qx(ij_u_we(k,1),ij_u_we(k,2)))*dx*dt
    end do
    do k=1,n_ij_v_sn
      qs_total=qs_total+abs(qsy(ij_v_sn(k,1),ij_v_sn(k,2)))*dx*dt
      q_total =q_total +abs( qy(ij_v_sn(k,1),ij_v_sn(k,2)))*dx*dt
    end do
    r_total=r_total+sum(r(:,:))*dx*dx*dt
    out_total=-out_total_ini+finf2_total+qs_total+q_total+eva1_total+eva2_total+(sum(h)+sum(hs*(lamda-pwc))+sum(pw*d))*dx*dx
    if(umax_disp<maxval(umax))then
      umax_disp=maxval(umax)
      maxloc_2d=maxloc(umax)
      i_umax=maxloc_2d(1)
      j_umax=maxloc_2d(2)
      h_umax_disp=h_umax(i_umax,j_umax)
    end if
    if(vmax_disp<maxval(vmax))then
      vmax_disp=maxval(vmax)
      maxloc_2d=maxloc(vmax)
      i_vmax=maxloc_2d(1)
      j_vmax=maxloc_2d(2)
      h_vmax_disp=h_vmax(i_vmax,j_vmax)
    end if
    ! <--- cal 2D
    !
    ! ---> cal 1d
    time2=time1
    itr2=1
    dtdiv=dt
    umax_1d=0.d0
    call cal_1d_uhb(ni_vect,i_vect,iend_1d,i_to,z_1d,h_1d,b_1d,n_1d,dtdiv,0.d0,dx_vect,u_1d,uo_1d,hlim_1d,umax_1d,q_1d)
    u_1d(i_outlet(:))=u_1d(i_from(i_outlet(:),1))
    flxd=umax_1d*dtdiv
    q_1d(i_outlet(:))=b_1d(i_outlet(:))*h_1d(i_outlet(:))*u_1d(i_outlet(:))
    q_in(i_scal(:))=q_1d(i_from(i_scal(:),1))+q_1d(i_from(i_scal(:),2))+q_1d(i_from(i_scal(:),3))
    q_in(i_inlet(1:ni_inlet))=0.d0 !qin
    q_out(:)=q_1d(i_node(:))
    dtmin_h1d=dtdiv
    ho_1d=h_1d
    call cal_1d_h(iend_1d,i_node,iend_1d,dtdiv,dx_scal,b_1d,h_1d,iend_1d,q_in,q_out,dtmin_h1d)
    h_1d=ho_1d
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ((dx*cfl) < flxd .or. dtmin_h1d < dt )then
      itr2=max(nint(dt/dtmin_h1d),ceiling(flxd/(dx*cfl)))
      if(itr2<=0)then
        write(*,*)'itr2',itr2,flxd,dx*cfl,nint(dt0/dtlim)
        read(*,*)
      end if
      dtdiv=dt/dble(itr2)
      if(dtdiv<dtmin)dtmin=dtdiv
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do while (itr2>0)
      dtmin_h1d=0.d0
      call cal_1d_h(iend_1d,i_node,iend_1d,dtdiv,dx_scal,b_1d,h_1d,iend_1d,q_in,q_out,dtmin_h1d)
      call cal_1d_uhb(ni_vect,i_vect,iend_1d,i_to,z_1d,h_1d,b_1d,n_1d,dtdiv,cfl,dx_vect,u_1d,uo_1d,hlim_1d,umax_1d,q_1d)
      u_1d(i_outlet(:))=u_1d(i_from(i_outlet(:),1))
      q_1d(i_outlet(:))=b_1d(i_outlet(:))*h_1d(i_outlet(:))*u_1d(i_outlet(:))
      q_in(i_scal(:))=qo_1d(i_from(i_scal(:),1))+qo_1d(i_from(i_scal(:),2))+qo_1d(i_from(i_scal(:),3))
      q_in(i_inlet(1:ni_inlet))=0.d0 !qin
      q_out(:)=qo_1d(i_node(:))
      uo_1d(:)=u_1d(:)
      ho_1d(:)=h_1d(:)
      qo_1d(:)=q_1d(:)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      q_1d_total=q_1d_total+sum(q_out(i_from(i_outlet(:),1)))*dtdiv
      h_1d(i_outlet(1:ni_outlet))=0.d0
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      time2=time2+dtdiv
      itr2=itr2-1
    end do ! div, do while
    ! <--- cal 1d
    !
    time1=time1+dt
    itr=itr-1
  end do ! do while
  time0=time
  dt=dt0
  ! <--- output
  if (mod(itsp, itsp_d_out) ==0) then
    ! ---> disp out
    write(*,*)
    write(*,'(i5,a,a,a,i5,a,a)')itsp/itsp_rain,' ',datetime(itsp/itsp_rain),' +', &
    int(mod(itsp,itsp_rain)*dt0),' (s) ---> ',datetime(itsp/itsp_rain+1)
    write(*,'(a,f10.4,a,f10.3)')'time(day) =',time/86400.d0,', time(hour)=',time/3600.d0
    write(*,'(a,i10,a,f10.3,a,f6.2,a,2i4)')&
          'time(s)   =',nint(time),', rain(mm/h)=',maxval(r(:,:))*1000.d0*3600.d0
    write(*,'(2(a,f6.2,a,2i4,a))')&
          ' qx_max=',maxval(abs(qx)) ,' @ i,j=',maxloc(abs(qx)) ,',', &
          ' qy_max=',maxval(abs(qy)) ,' @ i,j=',maxloc(abs(qy))
    write(*,'(2(a,f6.2,a,2i4,a))')&
          'qsx_max=',maxval(abs(qsx)),' @ i,j=',maxloc(abs(qsx)),',', &
          'qsy_max=',maxval(abs(qsy)),' @ i,j=',maxloc(abs(qsy))
    write(*,'(2(a,f6.2,a,2i4,a))')&
          '  h_max=',maxval(h)       ,' @ i,j=',maxloc(h)       ,',', &
          ' hs_max=',maxval(hs)      ,' @ i,j=',maxloc(hs)
    write(*,'(2(a,f10.3),a,2i5)')'umax=',umax_disp,', h@umax=',h_umax_disp,', i,j=',i_umax,j_umax
    write(*,'(2(a,f10.3),a,2i5)')'vmax=',vmax_disp,', h@vmax=',h_vmax_disp,', i,j=',i_vmax,j_vmax
    write(*,'(2(a,f10.5))')'dt_min=',dtmin,', dt_lim=',dtlim
    umax_disp=0.d0
    vmax_disp=0.d0
    dtmin=dt0
    write(*,'(4(a,e12.5))')&
          'rain=',r_total,',  pw=',sum(pw*d)*dx*dx,',  e1=',eva1_total,',  e2=',eva2_total
    write(*,'(4(a,e12.5))')&
          '  hs=',sum(hs*(lamda-pwc))*dx*dx,',   h=',(sum(h))*dx*dx,',  qs=',qs_total,',   q=',q_total
    write(*,'(3(a,e12.5))')' sum=',q_1d_total+out_total+sum(h_1d*b_1d*dx_scal), &
          ',h_1d=',sum(h_1d*b_1d*dx_scal),',q_1d=',q_1d_total
    !do i=1,ni_outlet
      !write(*,'(2(a,e12.5))')'q_1d=',q_1d_total, ',h_1d=',sum(h_1d*b_1d*dx_scal)
    !end do
    write(*,*)
    ! <--- disp out
    write(200,'(a,i10,9(a,E20.10e3),2(a,i10),6(a,E20.10e3))')&
          '  time(s)=',nint(time),',',r_total,',', &
    out_total,',',qs_total,',',q_total,',', &
    eva1_total,',',eva2_total,',',sum(hs*(lamda-pwc))*dx*dx,',',sum(h)*dx*dx,',', &
    sum(pw*d)*dx*dx,',',0,',',0,',', &
    finf2_total,',',q_1d_total,',',sum(h_1d*b_1d*dx_scal)
    write(201,'(a10,"",i10,*(e12.3,:,""))')'time(s)=',nint(time),vsl_out(:)
  end if
  !
  if (mod(itsp, itsp_f_out) ==0 .or. itsp==itsp_max) then
    write(fotime,'(i10.10)')nint(time)
    ! ---> 1d output
    fname=trim(adjustl(out_dir))//'res_1d_'//trim(adjustl(fotime))//'.txt'
    open(20,file=fname)
    do i=1,iend_1d
      write(20,'(2i5,i10,3f10.3,f15.3)')ij_2d1d(i,1),ij_2d1d(i,2),i,z_1d(i),h_1d(i),u_1d(i),q_1d(i)
    end do
    close(20)
    ! <--- 1d output
    cfmt='(*(f10.3))'
    fname=trim(adjustl(out_dir))//'res_2d_hs_'//trim(adjustl(fotime))//'.asc'
    call w_fasc(20,fname,cfmt,hs(:,:),iend,jend,xllcorner,yllcorner,dx,0.d0)
    !
    fname=trim(adjustl(out_dir))//'res_2d_h_'//trim(adjustl(fotime))//'.asc'
    call w_fasc(20,fname,cfmt,h(:,:),iend,jend,xllcorner,yllcorner,dx,0.d0)
    !
!    fname='output/pw_'//trim(adjustl(fotime))//'.asc'
!    call w_fasc(20,fname,cfmt,pw(:,:),iend,jend,xllcorner,yllcorner,dx,0.d0)
!    !
    fname=trim(adjustl(out_dir))//'res_2d_hmax.asc'
    call w_fasc(20,fname,cfmt,hmax(:,:),iend,jend,xllcorner,yllcorner,dx,0.d0)
!    !
    fname=trim(adjustl(out_dir))//'res_2d_hsmax.asc'
    call w_fasc(20,fname,cfmt,hsmax(:,:),iend,jend,xllcorner,yllcorner,dx,0.d0)
!    !
!    fname='output/qx_ave_'//trim(adjustl(fotime))//'.asc'
!    call w_fasc(20,fname,cfmt,qx_acc(:,:)/output_f_time,iend+1,jend+1, &
!    xllcorner-dx*0.5,yllcorner,dx,0.d0)
!    !
!    fname='output/qy_ave_'//trim(adjustl(fotime))//'.asc'
!    call w_fasc(20,fname,cfmt,qy_acc(:,:)/output_f_time,iend+1,jend+1, &
!    xllcorner,yllcorner-dx*0.5,dx,0.d0)
!    !
!    fname='output/qsx_ave_'//trim(adjustl(fotime))//'.asc'
!    call w_fasc(20,fname,cfmt,qsx_acc(:,:)/output_f_time,iend+1,jend+1, &
!    xllcorner-dx*0.5,yllcorner,dx,0.d0)
!    !
!    fname='output/qsy_ave_'//trim(adjustl(fotime))//'.asc'
!    call w_fasc(20,fname,cfmt,qsy_acc(:,:)/output_f_time,iend+1,jend+1, &
!    xllcorner,yllcorner-dx*0.5,dx,0.d0)
    !
    qx_acc(:,:)=0.0d0
    qy_acc(:,:)=0.0d0
    qsx_acc(:,:)=0.0d0
    qsy_acc(:,:)=0.0d0
    !
  end if
  ! <--- output
  !
end do
!
write(*,*)'--- nomal end ---'
!
stop
end
!
!
!
subroutine cal_q(n_1d,ij_1d, &
                 iend,jend,z,h,n,ioffset,joffset,hthr,dx,flux,uvmax,h_uvmax,dt,cfl)
use omp_lib
implicit none
real(8), parameter :: hm=5.0d0/3.0d0
integer :: ni,n_1d,ij_1d(1:n_1d,1:2),i,j,iend,jend,ip,jp,ioffset,joffset
real(8) :: z(iend,jend),h(iend,jend),wl(iend,jend),n(iend,jend), &
           flux(iend+1,jend+1),uvmax(iend+1,jend+1),h_uvmax(iend+1,jend+1), &
           grad,htmp,uv,dt,cfl,hthr,dx,cost,sint
!
!uvmax=0.d0
wl=z+h
!$omp parallel do private(ni,i,j,ip,jp,grad,cost,sint,htmp,uv)
do ni=1,n_1d
  i=ij_1d(ni,1)
  j=ij_1d(ni,2)
  ip=i+ioffset
  jp=j+joffset
  grad=(wl(i,j)-wl(ip,jp))/dx
  grad=sign(sin(atan(abs(grad))),grad)
  cost=cos(atan(grad))
  sint=sin(atan(grad))
  if(grad>0.0d0)then
    if(wl(i,j)<z(ip,jp) .or. h(i,j)<hthr)then
      flux(i,j)=0.0d0
    else
      htmp=wl(i,j)-max(z(i,j),z(ip,jp))
      htmp=htmp*cost
      flux(i,j)=sign(1.0d0/n(i,j)*abs(sint)**0.5d0*htmp**hm,-grad)
      uv=abs(flux(i,j))/htmp
    end if
  else
    if(wl(ip,jp)<z(i,j) .or. h(ip,jp)<hthr )then
      flux(i,j) = 0.0d0
    else
      htmp=wl(ip,jp)-max(z(i,j),z(ip,jp))
      htmp=htmp*cost
      flux(i,j)=sign(1.0d0/n(ip,jp)*abs(sint)**0.5d0*htmp**hm,-grad)
      uv=abs(flux(i,j))/htmp
    end if
  end if
  if(dx/dt*cfl < uv .and. 0.d0<cfl)then
    uv=sign(dx/dt*cfl,flux(i,j))
    flux(i,j)=uv*htmp
  end if
  if(uvmax(i,j)<uv)then
    uvmax(i,j)=uv
    h_uvmax(i,j)=htmp
  end if
end do
!$omp end parallel do
!call cal_q_bd(n_1d_bd,ij_1d_bd,iend,jend,z,h,n,dd,flux)
!
end subroutine cal_q
!
!
!
subroutine cal_q_bd(n_1d,ij_1d,iend,jend,z,h,n,dx,flux)
implicit none
real(8), parameter :: hm=5.0d0/3.0d0
integer :: ni,n_1d,ij_1d(n_1d,8),i,j,iend,jend,iofs1,jofs1,iofs2,jofs2,iofs3,jofs3
real(8) :: z(iend,jend),h(iend,jend),wl(iend,jend),n(iend,jend), &
           flux(iend+1,jend+1),grad(iend+1,jend+1),cost,sint,dx
!
! west  h(i  ,  j),(wl(i+1,  j)-wl(i  ,  j))/dx (0 , 0),( 1, 0),( 0, 0),flux<0
! east  h(i-1,  j),(wl(i-1,  j)-wl(i-2,  j))/dx (-1, 0),(-1, 0),(-2, 0),flux>0
! north h(i  ,j-1),(wl(i  ,j-1)-wl(i  ,j-2))/dy (0 ,-1),( 0,-1),( 0,-2),flux>0
! south h(i  ,  j),(wl(i  ,j+1)-wl(i  ,  j))/dy (0 , 0),( 0, 1),( 0, 0),flux<0
wl=z+h
do ni=1,n_1d
  i    =ij_1d(ni,1)
  j    =ij_1d(ni,2)
  iofs1=ij_1d(ni,3)
  jofs1=ij_1d(ni,4)
  iofs2=ij_1d(ni,5)
  jofs2=ij_1d(ni,6)
  iofs3=ij_1d(ni,7)
  jofs3=ij_1d(ni,8)
  grad(i,j)=(wl(i+iofs2,j+jofs2)-wl(i+iofs3,j+jofs3))/dx
  cost=cos(atan(grad(i,j)))
  sint=sin(atan(grad(i,j)))
  flux(i,j)=sign(1.0d0/n(i+iofs1,j+jofs1) *&
  abs(sint)**0.50d0*(h(i+iofs1,j+jofs1)*cost)**hm,-grad(i,j))
  if(iofs2>0 .or. jofs2>0)then
    if(flux(i,j)>0.0d0)flux(i,j)=0.0d0
  else if (iofs2<0 .or. jofs2<0)then
    if(flux(i,j)<0.0d0)flux(i,j)=0.0d0
  end if
end do
!
end subroutine cal_q_bd
!
!
!
subroutine cal_qs(n_1d,ij_1d, &
                  iend,jend,z,h,ioffset,joffset,kk,hthr,dx,flux,uvmax,h_uvmax)
use omp_lib
implicit none
integer :: ni,n_1d,ij_1d(1:n_1d,1:2),i,j,iend,jend,ip,jp,ioffset,joffset
real(8) :: z(iend,jend),h(iend,jend),wl(iend,jend),kk(iend,jend), &
           flux(iend+1,jend+1),uvmax(iend+1,jend+1),h_uvmax(iend+1,jend+1), &
           grad(iend+1,jend+1),cost,sint,hthr,dx,htmp,uv
!
!uvmax=0.d0
wl=z+h
!$omp parallel do private(ni,i,j,ip,jp,cost,sint,htmp,uv)
do ni=1,n_1d
  i=ij_1d(ni,1)
  j=ij_1d(ni,2)
  ip=i+ioffset
  jp=j+joffset
  grad(i,j)=(wl(i,j)-wl(ip,jp))/dx
  cost=cos(atan(grad(i,j)))
  sint=sin(atan(grad(i,j)))
  if ( grad(i,j) > 0.0d0 ) then
    htmp=wl(i,j) -max(z(i,j),z(ip,jp))
    htmp=htmp*cost
    flux(i,j)=sign(abs(sint)*(wl(i,j)-max(z(i,j),z(ip,jp)))*kk(i,j),-grad(i,j))
    if(htmp>0.d0)then
      uv=abs(flux(i,j))/htmp
    else
      uv=0.d0
    end if
    if ( wl(i,j) < z(ip,jp) .or. h(i,j) < hthr ) then
      flux(i,j) = 0.0d0
      uv = 0.d0
    end if
  else
    htmp=wl(ip,jp)-max(z(i,j),z(ip,jp))
    htmp=htmp*cost
    flux(i,j)=sign(abs(sint)*(wl(ip,jp)-max(z(i,j),z(ip,jp)))*kk(ip,jp), &
    -grad(i,j))
    if(htmp>0.d0)then
      uv=abs(flux(i,j))/htmp
    else
      uv=0.d0
    end if
    if ( wl(ip,jp) < z(i,j) .or. h(ip,jp) < hthr )then
      flux(i,j) = 0.0d0
      uv = 0.d0
    end if
  end if
  if(uvmax(i,j)<uv)then
    uvmax(i,j)=uv
    h_uvmax(i,j)=htmp
  end if
end do
!$omp end parallel do
!call cal_qs_bd(n_1d_bd,ij_1d_bd,iend,jend,z,h,kk,dd,flux)
!
end subroutine cal_qs
!
!
!
subroutine cal_qs_bd(n_1d,ij_1d,iend,jend,z,h,kk,dx,flux)
implicit none
integer :: ni,n_1d,ij_1d(n_1d,8),i,j,iend,jend,iofs1,jofs1,iofs2,jofs2,iofs3,jofs3
real(8) :: z(iend,jend),h(iend,jend),wl(iend,jend),kk(iend,jend), &
           flux(iend+1,jend+1),grad(iend+1,jend+1),cost,sint,dx
!
! west  h(i  ,  j),(wl(i+1,  j)-wl(i  ,  j))/dx (0 , 0),( 1, 0),( 0, 0),flux<0
! east  h(i-1,  j),(wl(i-1,  j)-wl(i-2,  j))/dx (-1, 0),(-1, 0),(-2, 0),flux>0
! north h(i  ,j-1),(wl(i  ,j-1)-wl(i  ,j-2))/dy (0 ,-1),( 0,-1),( 0,-2),flux>0
! south h(i  ,  j),(wl(i  ,j+1)-wl(i  ,  j))/dy (0 , 0),( 0, 1),( 0, 0),flux<0
wl=z+h
do ni=1,n_1d
  i    =ij_1d(ni,1)
  j    =ij_1d(ni,2)
  iofs1=ij_1d(ni,3)
  jofs1=ij_1d(ni,4)
  iofs2=ij_1d(ni,5)
  jofs2=ij_1d(ni,6)
  iofs3=ij_1d(ni,7)
  jofs3=ij_1d(ni,8)
  grad(i,j)=(wl(i+iofs2,j+jofs2)-wl(i+iofs3,j+jofs3))/dx
  cost=cos(atan(grad(i,j)))
  sint=sin(atan(grad(i,j)))
  flux(i,j)=sign(abs(sint)*(h(i+iofs1,j+jofs1)*cost)*kk(i+iofs1,j+jofs1), &
  -grad(i,j))
  if(iofs2>0 .or. jofs2>0)then
    if(flux(i,j)>0.0d0)flux(i,j)=0.0d0
  else if (iofs2<0 .or. jofs2<0)then
    if(flux(i,j)<0.0d0)flux(i,j)=0.0d0
  end if
end do
!
end subroutine cal_qs_bd
!
!
!
subroutine cal_divha(n_1d,ij_1d,iend,jend,lamda,pwc,ha,h,hs,pw,d,finf2,dt,qx,qy,qsx,qsy)
use omp_lib
implicit none
integer :: ni,n_1d,ij_1d(n_1d,2),i,j,iend,jend
real(8) :: lamda(iend,jend),pwc(iend,jend),ho, &
           ha(iend,jend),h(iend,jend),hs(iend,jend),pw(iend,jend),d(iend,jend), &
           finf2(iend,jend),dt,qx(iend+1,jend+1),qy(iend+1,jend+1), &
           qsx(iend+1,jend+1),qsy(iend+1,jend+1)
!
!$omp parallel do private(ni,i,j,ho)
do ni=1,n_1d
  i=ij_1d(ni,1)
  j=ij_1d(ni,2)
  ho=h(i,j)
  if(ha(i,j)>lamda(i,j)*d(i,j)+finf2(i,j)*dt)then
    h(i,j)=ha(i,j)-(lamda(i,j)*d(i,j)+finf2(i,j)*dt)
    hs(i,j)=d(i,j)
    pw(i,j)=pwc(i,j)
  else if(ha(i,j)>pwc(i,j)*d(i,j)+finf2(i,j)*dt)then
    h(i,j)=0.d0
    hs(i,j)=(ha(i,j)-(pwc(i,j)*d(i,j)+finf2(i,j)*dt))/(lamda(i,j)-pwc(i,j))
    pw(i,j)=pwc(i,j)
  else if(ha(i,j)>pwc(i,j)*d(i,j))then
    h(i,j)=0.d0
    hs(i,j)=0.d0
    pw(i,j)=pwc(i,j)
    finf2(i,j)=(ha(i,j)-pwc(i,j)*d(i,j))/dt
  else if(ha(i,j)>=0.d0)then
    h(i,j)=0.d0
    hs(i,j)=0.d0
    pw(i,j)=ha(i,j)/d(i,j)
    finf2(i,j)=0.d0
  else
    write(*,'(a)')'ha < 0, @cal_divha -> stop'
    write(*,'(2i5,2f10.3)')i,j,ha(i,j),ho
    write(*,'(5f10.3)')qx(i+1,j),qx(i,j),qy(i,j+1),qy(i,j), &
                  &  (-(qx(i+1,j)-qx(i,j))-(qy(i,j+1)-qy(i,j)))*dt
    write(*,'(5f10.3)')qsx(i+1,j),qsx(i,j),qsy(i,j+1),qsy(i,j), &
                  &  (-(qsx(i+1,j)-qsx(i,j))-(qsy(i,j+1)-qsy(i,j)))*dt
    if(ha(i,j)<-0.01d0)then
    read(*,*)
    end if
    h(i,j)=0.d0
    ha(i,j)=0.d0
    !read(*,*)
    !stop
  end if
end do
!$omp end parallel do
!
end subroutine cal_divha
!
!
!
subroutine cal_eva(n_1d,ij_1d,iend,jend,r,hs,pw,lamda,pwc,d,dt,eva1,eva2)
use omp_lib
implicit none
integer :: ni,n_1d,ij_1d(n_1d,2),i,j,iend,jend
real(8) :: lamda(iend,jend),pwc(iend,jend),dt,d(iend,jend), &
           r(iend,jend),hs(iend,jend),pw(iend,jend),eva1(iend,jend),eva2(iend,jend)
!
!$omp parallel do private(ni,i,j)
do ni=1,n_1d
  i=ij_1d(ni,1)
  j=ij_1d(ni,2)
  if(r(i,j)>0.d0)then
    eva1(i,j)=0.d0
    eva2(i,j)=0.d0
  else if(hs(i,j)>0.d0)then
    eva1(i,j)=min(eva1(i,j),(lamda(i,j)-pwc(i,j))*hs(i,j)/dt)
    hs(i,j)=hs(i,j)-eva1(i,j)*dt/(lamda(i,j)-pwc(i,j))
    eva2(i,j)=0.d0
  else if(pw(i,j)>=0.d0)then
    eva1(i,j)=0.d0
    eva2(i,j)=min(eva2(i,j),pw(i,j)*d(i,j)/dt)
    pw(i,j)=pw(i,j)-eva2(i,j)*dt/d(i,j)
  end if
end do
!$omp end parallel do
!
end subroutine cal_eva
!
!
!
subroutine w_fasc(fon,fname,cfmt,a,iend,jend,xllcorner,yllcorner,cellsize,nv)
implicit none
integer :: j,iend,jend,fon
real(8) :: a(iend,jend),xllcorner,yllcorner,cellsize,nv
character(len=100) :: cfmt,fname
!
open(fon,file=trim(adjustl(fname)))
write(fon,'(a,i5)')'ncols',iend
write(fon,'(a,i5)')'nrows',jend
write(fon,'(a,f15.3)')'xllcorner',xllcorner
write(fon,'(a,f15.3)')'yllcorner',yllcorner
write(fon,'(a,f10.3)')'cellsize ',cellsize
write(fon,'(a,f10.3)')'NODATA_value',nv
do j=jend,1,-1
  write(fon,trim(cfmt))a(1:iend,j)
end do
close(fon)
!
end subroutine w_fasc
!
!
!
subroutine w_iasc(fon,fname,ia,iend,jend,xllcorner,yllcorner,cellsize,nv)
implicit none
integer :: j,iend,jend,fon,nv,ia(iend,jend)
real(8) :: xllcorner,yllcorner,cellsize
character(len=100) :: fname
!
open(fon,file=trim(adjustl(fname)))
write(fon,'(a,i5)')'ncols',iend
write(fon,'(a,i5)')'nrows',jend
write(fon,'(a,f15.3)')'xllcorner',xllcorner
write(fon,'(a,f15.3)')'yllcorner',yllcorner
write(fon,'(a,f10.3)')'cellsize ',cellsize
write(fon,'(a,i5)')'NODATA_value',nv
do j=jend,1,-1
  write(fon,'(*(i2))')ia(1:iend,j)
end do
close(fon)
!
end subroutine w_iasc
!
!
!
subroutine r_iasc(fon,fname,ia,iend,jend,xllcorner,yllcorner,cellsize)
implicit none
integer :: j,iend,jend,fon,ncols,nrows,ia(iend,jend)
real(8) :: xllcorner,yllcorner,cellsize
character(len=100) :: ctmp,fname
!
open(fon,file=trim(adjustl(fname)))
read(fon,*)ctmp,ncols
read(fon,*)ctmp,nrows
if (ncols /= iend .or. nrows /= jend) then
  write(*,*)'Check ncols /= iend or nrows /= jend'
  stop
end if
read(fon,*)ctmp,xllcorner
read(fon,*)ctmp,yllcorner
read(fon,*)ctmp,cellsize
read(fon,*)
do j=jend,1,-1
  read(fon,*)ia(1:iend,j)
end do
close(fon)
!
end subroutine r_iasc
!
!
!
subroutine r_fasc(fon,fname,a,iend,jend,xllcorner,yllcorner,cellsize)
implicit none
integer :: j,iend,jend,fon,ncols,nrows
real(8) :: xllcorner,yllcorner,cellsize, a(iend,jend)
character(len=100) :: ctmp,fname
!
open(fon,file=trim(adjustl(fname)))
read(fon,*)ctmp,ncols
read(fon,*)ctmp,nrows
if (ncols /= iend .or. nrows /= jend) then
  write(*,*)'Check ncols /= iend or nrows /= jend'
  stop
end if
read(fon,*)ctmp,xllcorner
read(fon,*)ctmp,yllcorner
read(fon,*)ctmp,cellsize
read(fon,*)
do j=jend,1,-1
  read(fon,*)a(1:iend,j)
end do
close(fon)
!
end subroutine r_fasc
!
!
!
subroutine point(fon,fname,a,iend,jend,xllcorner,yllcorner,cellsize)
implicit none
integer :: i,j,iend,jend,fon,a(iend,jend),iflg
real(8) :: xllcorner,yllcorner,cellsize,x,y
character(len=100) :: fname
!
open(fon,file=trim(adjustl(fname)))
write(fon,'(a)')'x,y,i,j,flag'
do i=1,iend
  do j=1,jend
    x=xllcorner+cellsize*float(i-1)
    y=yllcorner+cellsize*float(j-1)
    iflg=int(a(i,j))
    write(fon,'(2(f15.3,a),3(i5,a))')x,',',y,',',i,',',j,',',iflg
  end do
end do
close(fon)
!
end subroutine point
!
!
!
subroutine point2(fon,fname,xllcorner,yllcorner,cellsize,n_1d,i_1d,j_1d)
implicit none
integer :: ni,n_1d,i_1d(1:n_1d),j_1d(1:n_1d),i,j,fon
real(8) :: xllcorner,yllcorner,cellsize
character(len=100) :: fname
!
open(fon,file=trim(adjustl(fname)))
write(fon,'(a)')'x,y,i,j'
do ni=1,n_1d
  i=i_1d(ni)
  j=j_1d(ni)
  write(fon,'(2(f15.3,a),2(i5,a))')&
  xllcorner+cellsize*float(i-1),',',yllcorner+cellsize*float(j-1),',',i,',',j
end do
close(fon)
!
end subroutine point2
!
!
!
subroutine cal_dq(n_1d,ij_1d,iend,jend,dx,h,qx,qy,dq,dt)
use omp_lib
implicit none
integer :: ni,n_1d,ij_1d(n_1d,2),i,j,iend,jend
real(8) :: dx,dt,h(iend,jend),dq(iend,jend), &
           qx(iend+1,jend+1),qy(iend+1,jend+1)
!
!$omp parallel do private(ni,i,j)
do ni=1,n_1d
  i=ij_1d(ni,1)
  j=ij_1d(ni,2)
  dq(i,j)=-((qx(i+1,j)-qx(i,j))/dx+(qy(i,j+1)-qy(i,j))/dx)
  if(h(i,j) + dq(i,j)*dt < 0.d0 .and. 0.d0 < dt)then
    dt=min(dt,h(i,j)/dq(i,j))
  end if
end do
!$omp end parallel do
!
end subroutine cal_dq
!
!
!
subroutine cal_1d_uhb(ni_vect,i_vect,iend,i_to,z,h,b,n,dt,cfl,dx,u,uo,hlim,umax,q)
use omp_lib
implicit none
real(8), parameter :: PI=acos(-1.d0),deg2rad=PI/180.d0,rad2deg=180.d0/PI, &
                      hm=5.d0/3.d0,rm=4.d0/3.d0
integer :: i,iend,j,ni_vect,i_vect(ni_vect),i_to(iend)
real(8) :: hlim,u_tmp,h1,r1,b1,umax,cost,sint,dt,cfl,grad(iend),dx(iend), &
           z(iend),h(iend),u(iend),uo(iend),n(iend),b(iend),q(0:iend)
!
!        from  i    to
! FX      q    q    q
!         ->   ->   ->
! CV | z  | z  | z  |
!    |from| i  | to |
!umax=0.d0
!$omp parallel do private(j,i,cost,sint,h1,r1,b1,u_tmp) reduction(max: umax)
do j=1,ni_vect
  i=i_vect(j)
  if(h(i)<hlim .and. h(i_to(i))<hlim)then
    u(i)=0.d0
    q(i)=0.d0
    cycle
  end if
  grad(i)=(z(i_to(i))+h(i_to(i))-z(i)-h(i))/dx(i)
  cost=cos(atan(grad(i)))
  sint=sin(atan(grad(i)))
  !if(uo(i)>=0.d0)then ! Be careful when uo=0
  if(uo(i)>0.d0 .or. (uo(i)==0.d0 .and. grad(i)<=0.d0 .and. h(i)>=hlim))then
    h1=z(i)+h(i)-max(z(i),z(i_to(i)))
    h1=h1*cost
    r1=h1*b(i)/(2.d0*h1+b(i))
    b1=b(i)
    if(h1<hlim)then
      u(i)=0.d0
      q(i)=0.d0
      cycle
    end if
  else if(uo(i)<0.d0 .or. (uo(i)==0.d0 .and. grad(i)>0.d0 .and. h(i_to(i))>=hlim))then
    h1=z(i_to(i))+h(i_to(i))-max(z(i),z(i_to(i)))
    h1=h1*cost
    r1=h1*b(i_to(i))/(2.d0*h1+b(i_to(i)))
    b1=b(i_to(i))
    if(h1<hlim)then
      u(i)=0.d0
      q(i)=0.d0
      cycle
    end if
  else
    u(i)=0.d0
    q(i)=0.d0
    cycle
  end if
  !u_tmp=sign(1.d0/rn*h1**(2.d0/3.d0)*abs(grad(i))**0.5d0,-grad(i))
  u_tmp=sign(1.d0/n(i)*r1**(2.d0/3.d0)*abs(sint)**0.5d0,-grad(i))
  u(i)=u_tmp
  if(uo(i)*u(i) < 0.d0)u(i)=0.d0
  if(dx(i)/dt*cfl<abs(u(i)) .and. 0.d0<cfl)then
 !   write(*,*)'@u_1d cfl',u(i),dx(i)/dt*cfl
 !   write(*,*)'@u_1d cfl',dt,dx(i)*cfl/abs(u(i)),dt-dx(i)*cfl/abs(u(i))
    !read(*,*)
    u(i)=sign(dx(i)/dt*cfl,u(i))
  end if
  q(i)=b1*h1*u(i) ! A*V
  if(abs(u(i))>umax)umax=abs(u(i))
end do
!$omp end parallel do
!
end subroutine cal_1d_uhb
!
!
!
subroutine cal_1d_h(ni_scal,i_scal,iend,dt,dx,b,h,iqnum,q_in,q_out,dtmin_h1d)
use omp_lib
implicit none
real(8),parameter :: dtmin0=0.00001d0
integer :: i,iend,j,ni_scal,i_scal(ni_scal),iqnum
real(8) :: dt,dt2,dtmin_h1d,b(iend),h(iend),ho,dx(iend),q_in(iqnum),q_out(iqnum)
!
!        from  i    to
! FX      q    q    q
!         ->   ->   ->
! CV | z  | z  | z  |
!    |from| i  | to |
dt2=dt
!$omp parallel do private(j,i,ho) reduction(min: dt2,dtmin_h1d)
do j=1,ni_scal
  i=i_scal(j)
  ho=h(i)
  h(i)=ho-(q_out(i)-q_in(i))/dx(i)/b(i)*dt
  if(h(i)<0.d0)then
    if(dtmin_h1d<=0.d0)then
      write(*,'(a,4f15.7)')'h<0 @h_1d (dt,h,qo,qi)=',dt,h(i),q_out(i),q_in(i)
      !h(i)=0.d0
      !q_out(i)=0.d0
      read(*,*)
    end if
    !write(*,'(a,2i5,5f10.5)')'@1dh h1d < 0 ',i,j,ho,h(i),q_out(i),q_in(i)
    !write(*,'(a,f10.5)')'dt=',dt
    !read(*,*)
    !h(i)=0.d0
    !q_out(i)=0.d0
    if(0.d0 < (q_out(i)-q_in(i)))then
      dt2=min(dt2,ho/(q_out(i)-q_in(i))*dx(i)*b(i))
    else
    !  write(*,*)'qout-qin=0'
    !  read(*,*)
      !dt2=dt
    end if
    !write(*,'(a,3f15.7)')'@h_1d (dt2,h2,dt)=',dt2,h(i),dt
    if(dt2<0.00001d0)then
      write(*,'(a,3f15.7)')'dt<lim@h_1d (dt2,h2,dt)=',dt2,h(i),dt
      write(*,*)'stop'
      read(*,*)
    end if
  dtmin_h1d=dt2
  end if
end do
!$omp end parallel do
!
end subroutine cal_1d_h
!
!
!
subroutine cal_2d1d(iend,jend,ij_2d1d,iend_1d,dx,dtdiv,z,h,z_1d,h_1d,b_1d,dx_scal,vsl_out)
!https://www.mlit.go.jp/river/shishin_guideline/gijutsu/gijutsukijunn/chousa/pdf/07.pdf
use omp_lib
implicit none
real(8), parameter :: u1=2.d0/3.d0**(3.d0/2.d0),u2=0.35d0,u3=0.91d0,g=9.8d0
integer :: i,j,k,iend,jend,iend_1d,ij_2d1d(iend_1d,2)
real(8) :: dx,dtdiv,wl2d,wl1d,h1,h2,qsr,vsr,z(iend,jend),h(iend,jend), &
           z_1d(iend_1d),h_1d(iend_1d),b_1d(iend_1d),dx_scal(iend_1d),vsl_out(iend_1d)
!
!$omp parallel do private(k,i,j,wl2d,wl1d,h1,h2,qsr,vsr)
do k=1,iend_1d
  i=ij_2d1d(k,1)
  j=ij_2d1d(k,2)
  wl2d=z(i,j)+h(i,j)
  wl1d=z_1d(k)+h_1d(k)
  if(wl1d <= wl2d)then ! s->r from surface
    h1=h(i,j)
    h2=max(0.d0,wl1d-z(i,j))
    !if(z(i,j) <= wl1d)then
      if(h2<2.d0/3.d0*h1)then
        qsr=u2*h1*(2.d0*g*h1)**0.5d0
      else
        qsr=u3*h2*(2.d0*g*(h1-h2))**0.5d0
      end if
    !else                 ! fall
    !  qsr=u1*h1*(g*h1)**0.5d0
    !end if
  else ! r->s from surface
    h1=wl1d-z(i,j)
    h2=h(i,j)
    if(h2/h1<2.d0/3.d0)then
      qsr= -u2*h1*(2.d0*g*h1)**0.5d0
    else
      qsr= -u3*h2*(2.d0*g*(h1-h2))**0.5d0
    end if
  end if
  vsr=qsr*dx_scal(k)*dtdiv*2.d0
  if(vsr>=0.d0)then
    vsr=min(vsr,h(i,j)*(dx*dx))
  else
    vsr=max(vsr,-(wl1d-z(i,j))*(b_1d(k)*dx_scal(k)))
  end if
  vsl_out(k)=vsl_out(k)+vsr
  h_1d(k)=h_1d(k)+vsr/(b_1d(k)*dx_scal(k))
  h(i,j)=h(i,j)-vsr/(dx*dx)
  if(h(i,j)<0.d0)then
    !write(*,*)'@2d1d h2d<0',h(i,j),i,j
    !read(*,*)
    h(i,j)=0.d0
  end if
  if(h_1d(k)<0.d0)then
    !write(*,*)'@2d1d h1d<0',h_1d(k),k
    !read(*,*)
    h_1d(k)=0.d0
  end if
end do
!$omp end parallel do
!
end subroutine cal_2d1d
!
!
!
subroutine cal_grad(iend,jend,dx,z,grdave)
implicit none
integer :: i,iend,j,jend
real(8) :: dx,z(iend,jend),dzdx(iend,jend),dzdy(iend,jend),grdave(iend,jend)
!
dzdx(2:iend-1,:)=(z(3:iend,:)-z(1:iend-2,:))/(2.d0*dx)
dzdx(1,:)=(z(2,:)-z(1,:))/dx
dzdx(iend,:)=(z(iend,:)-z(iend-1,:))/dx
!
dzdy(:,2:jend-1)=(z(:,3:jend)-z(:,1:jend-2))/(2.d0*dx)
dzdy(:,1)=(z(:,2)-z(:,1))/dx
dzdy(:,jend)=(z(:,jend)-z(:,jend-1))/dx
!
grdave(:,:)=(dzdx**2.d0+dzdy**2.d0)**0.5d0
!
end subroutine cal_grad
