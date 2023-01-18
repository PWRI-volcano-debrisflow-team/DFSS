! Program: DR_ver_1.0.f90
! Copyright (2021) by Public Works Research Institute (P.W.R.I.)
! License: CC-BY-SA
!-------------------------------------------------------------------------------
!     1    2    3      end  end+1
! FX  q    q    q       q   q
!     ->   ->  ->   ->  ->  ->
! CV  | z  | z  | z |...| z |
!     | 1  | 2  | 3 |...|end|
!-------------------------------------------------------------------------------
implicit none
real(8), parameter :: PI=acos(-1.d0),deg2rad=PI/180.d0,rad2deg=180.d0/PI, &
                      g=9.8d0,u2=0.35d0,dtlim=0.0001d0,cfl=0.01d0
real(8), allocatable :: z(:,:),zo(:,:),zini(:,:),zs(:,:),h(:,:),ho(:,:),d(:,:), &
                        hs(:,:),hso(:,:),coh(:,:),sf(:,:),hmax(:,:),hsmax(:,:), &
                        cc(:,:),cco(:,:),cf(:,:),cfo(:,:),e(:,:), tant_e(:,:),&
                        grad_z(:,:),grad_w(:,:),rho(:,:),rhom(:,:),fdir(:,:)
real(8), allocatable :: u(:,:),uo(:,:),v(:,:),vo(:,:),ua(:,:),va(:,:),uv(:,:), &
                        uh(:,:),vh(:,:),qx(:,:),qy(:,:),qsx(:,:),qsy(:,:), &
                        qccx(:,:),qccy(:,:),qcfx(:,:),qcfy(:,:), &
                        q_a(:,:,:),qcc_a(:,:,:),qcf_a(:,:,:)
real(8), allocatable :: qin0(:),cin0(:)
integer, allocatable :: iacc(:,:),ibasin(:,:),idir(:,:),ifld(:,:)
real(8) :: dx,dy,dt,qdt,dout_time,fout_time,xllcorner,yllcorner,cellsize, &
           dm,csta,cstad,pf,pc,rhow,sig,tanp,coh_tmp,ks,rn,hlim
real(8) :: time,time0,time1,end_time,dt0,dtmin, &
           x,y,flxd,umax,vmax,h_umax,h_vmax,tmp,usmax,vsmax, &
           qin,cin,q_out,qcc_out,qcf_out
integer :: iend,jend,itsp_max,itsp_dout,itsp_fout,qnum,itsp_q
integer :: i,j,k,itsp,itr,i_umax,j_umax,i_vmax,j_vmax
integer :: ijofs_uv(8,0:4,2)= &    ! see subroutine cal_u
reshape((/0, 0,-1,-1, 0, 1, 1, 0,& ! @u case0 i=1:4,j=5:8
          0,-1, 0, 0, 0, 0, 0,-1,& ! @u case1 i=1:4,j=5:8
          0,-1, 0, 0, 0, 0, 1, 0,& ! @u case2 i=1:4,j=5:8
          1, 0, 0, 0, 0, 0, 0,-1,& ! @u case3 i=1:4,j=5:8
          1, 0, 0, 0, 0, 0, 1, 0,& ! @u case4 i=1:4,j=5:8
          0, 1, 1, 0, 0, 0,-1,-1,& ! @v case0 i=1:4,j=5:8
          0, 0, 0,-1, 0,-1, 0, 0,& ! @v case1 i=1:4,j=5:8
          0, 0, 1, 0, 0,-1, 0, 0,& ! @v case2 i=1:4,j=5:8
          0, 0, 0,-1, 1, 0, 0, 0,& ! @v case3 i=1:4,j=5:8
          0, 0, 1, 0, 1, 0, 0, 0 & ! @v case4 i=1:4,j=5:8
          /),(/8,5,2/))            ! 0:4 ---> 5
integer, allocatable :: ij_cv(:,:),ij_u(:,:),ij_v(:,:),ij_u_we(:,:),ij_v_sn(:,:), &
                        ij_w(:,:),ij_e(:,:),ij_s(:,:),ij_n(:,:)
integer :: n_ij_cv,n_ij_u,n_ij_v,n_ij_u_we,n_ij_v_sn, &
           n_ij_w,n_ij_e,n_ij_s,n_ij_n
! ---> 1d
integer,allocatable :: ij_1d2d(:,:,:),ivin(:),i_1d(:,:)
real(8),allocatable :: xy_1d2d(:,:),va_in0(:,:),vcc_in0(:,:),vcf_in0(:,:), &
                       va_in(:),vcc_in(:),vcf_in(:)
integer :: iend_1d,itmp,n_ivin,itr_vin,io,n_ivin0,itmp1,itmp2
real(8) :: va_acc,vcc_acc,vcf_acc
! <--- 1d
real(8), allocatable :: rseries(:,:,:),r(:,:),rmap(:,:,:),x_rpnt(:),y_rpnt(:),r_rpnt(:),d_idw(:)
real(8) :: rdt,rain_tmp
integer :: l,rnum,itsp_r,n_rpnt,minl(1)
character(len=100) :: rname
!
integer :: nsci,nscj
integer, allocatable :: iffx(:,:),iffy(:,:)
real(8),allocatable :: mvals(:)
character :: cline*1000
!
character(len=100) :: fname,ctmp,cfmt,out_dir
!
open(1,file='DF_input.txt')
read(1,*) ! grid
read(1,*)iend
read(1,*)jend
read(1,*)dx
read(1,*)dy
!
allocate(z(iend,jend),zo(iend,jend),zini(iend,jend),zs(iend,jend),h(iend,jend),ho(iend,jend),d(iend,jend),&
         hs(iend,jend),hso(iend,jend),coh(iend,jend),sf(iend,jend),hmax(iend,jend),hsmax(iend,jend), &
         cc(iend,jend),cco(iend,jend),cf(iend,jend),cfo(iend,jend),e(iend,jend),tant_e(iend,jend), &
         grad_z(iend,jend),grad_w(iend,jend),rho(iend,jend),rhom(iend,jend),fdir(iend,jend))
allocate(u(iend,jend),uo(iend,jend),v(iend,jend),vo(iend,jend),ua(iend,jend),va(iend,jend),uv(iend,jend), &
         uh(iend,jend),vh(iend,jend),qx(iend,jend),qy(iend,jend),qsx(iend,jend),qsy(iend,jend), &
         qccx(iend,jend),qccy(iend,jend),qcfx(iend,jend),qcfy(iend,jend), &
         q_a(iend,jend,2),qcc_a(iend,jend,2),qcf_a(iend,jend,2))
allocate(qin0(0:qnum),cin0(0:qnum))
allocate(iacc(iend,jend),ibasin(iend,jend),idir(iend,jend),ifld(iend,jend))
!
read(1,*)
read(1,*)dt
read(1,*)dout_time
read(1,*)fout_time
read(1,*)
read(1,*)rname
read(1,*)rnum
read(1,*)rdt
read(1,*)
read(1,*)qnum
read(1,*)qdt
read(1,*)
read(1,*)fname ! elevation
call r_fasc(10,fname,z(:,:),iend,jend,xllcorner,yllcorner,cellsize)
!
read(1,*)fname ! targetArea
call r_iasc(10,fname,ibasin(:,:),iend,jend,xllcorner,yllcorner,cellsize)
!
ifld(:,:)=0
!fname='input/fld.asc'
!call r_iasc(10,fname,ifld(:,:),iend,jend,xllcorner,yllcorner,cellsize)
!open(11,file='ifld_in.txt')
!read(11,*)nifld
!do k=1,nifld
!  read(11,*)i,j
!  ifld(i,j)=1
!end do
!
read(1,*)fname ! flowDir
call r_iasc(10,fname,idir(:,:),iend,jend,xllcorner,yllcorner,cellsize)
fdir(:,:)=idir(:,:)*45*deg2rad
!
read(1,*)fname ! soil depth(m)
call r_fasc(10,fname,d(:,:),iend,jend,xllcorner,yllcorner,cellsize)
zs(:,:)=z(:,:)-d(:,:)
zini(:,:)=z(:,:)
!
read(1,*)fname ! initial hs(m)
call r_fasc(10,fname,hs,iend,jend,xllcorner,yllcorner,cellsize)
hso(:,:)=hs(:,:)
!
read(1,*)
read(1,*)dm
read(1,*)csta
read(1,*)pf
read(1,*)rhow
read(1,*)sig
read(1,*)tanp ! deg
tanp=tan(tanp*deg2rad) ! deg ---> rad
read(1,*)coh_tmp
read(1,*)ks ! (mm/h)
ks=ks/100.0d0 ! (mm/h) -> (m/s)
read(1,*)rn
read(1,*)hlim
!
cstad=csta
pc=1.d0-pf
end_time=dble(qnum-1)*qdt
itsp_max=nint(end_time/dt)
itsp_dout=nint(dout_time/dt)
itsp_fout=nint(fout_time/dt)
itsp_q=nint(qdt/dt)
itsp_r=nint(rdt/dt)
out_dir='output_DF/'
write(*,*)trim(adjustl(out_dir))
write(*,*)''
!
read(1,*)
read(1,*)nsci
allocate(iffx(nsci,3))
write(cline,'(a10)')'time'
do i=1,nsci
  read(1,*)iffx(i,1),iffx(i,2),iffx(i,3) ! i,j1,j2
  write(ctmp,'(a,i0,a,i0,a,i0,a)')'i(',iffx(i,1),')&j(',iffx(i,2),':',iffx(i,3),'),'
  ctmp=adjustr(ctmp)
  cline=trim(cline)//ctmp(56:100)
end do
!
read(1,*)nscj
allocate(iffy(nscj,3))
do i=1,nscj
  read(1,*)iffy(i,1),iffy(i,2),iffy(i,3)
  write(ctmp,'(a,i0,a,i0,a,i0,a)')'i(',iffy(i,2),':',iffy(i,3),')&j(',iffy(i,1),'),'
  ctmp=adjustr(ctmp)
  cline=trim(cline)//ctmp(56:100)
end do
open(30,file='flux_sed.txt')
write(30,'(a)')cline
!
allocate(mvals((nsci+nscj)*3))
!
! ---> set cv & flux
open(1,file='input/controlVolume_centerPoint_inDF.txt')
read(1,*)n_ij_cv
allocate(ij_cv(n_ij_cv,2))
read(1,*)
do k=1,n_ij_cv
  read(1,*)x,y,ij_cv(k,1:2)
end do
close(1)
!
open(1,file='input/fluxPoint_x_inDF.txt')
read(1,*)n_ij_u
allocate(ij_u(n_ij_u,2))
read(1,*)
do k=1,n_ij_u
  read(1,*)x,y,ij_u(k,1:2)
end do
close(1)
!
open(1,file='input/fluxPoint_y_inDF.txt')
read(1,*)n_ij_v
allocate(ij_v(n_ij_v,2))
read(1,*)
do k=1,n_ij_v
  read(1,*)x,y,ij_v(k,1:2)
end do
close(1)
!
open(1,file='input/fluxPoint_xBoudary_inDF.txt')
read(1,*)n_ij_u_we
allocate(ij_u_we(n_ij_u_we,2))
read(1,*)
do k=1,n_ij_u_we
  read(1,*)x,y,ij_u_we(k,1:2)
end do
close(1)
!
open(1,file='input/fluxPoint_yBoudary_inDF.txt')
read(1,*)n_ij_v_sn
allocate(ij_v_sn(n_ij_v_sn,2))
read(1,*)
do k=1,n_ij_v_sn
  read(1,*)x,y,ij_v_sn(k,:)
end do
close(1)
!
open(1,file='input/fluxPoint_wBoudary_inDF.txt')
read(1,*)n_ij_w
allocate(ij_w(n_ij_w,2))
read(1,*)
do k=1,n_ij_w
  read(1,*)x,y,ij_w(k,:)
end do
close(1)
!
open(1,file='input/fluxPoint_eBoudary_inDF.txt')
read(1,*)n_ij_e
allocate(ij_e(n_ij_e,2))
read(1,*)
do k=1,n_ij_e
  read(1,*)x,y,ij_e(k,:)
end do
close(1)
!
open(1,file='input/fluxPoint_sBoudary_inDF.txt')
read(1,*)n_ij_s
allocate(ij_s(n_ij_s,2))
read(1,*)
do k=1,n_ij_s
  read(1,*)x,y,ij_s(k,:)
end do
close(1)
!
open(1,file='input/fluxPoint_nBoudary_inDF.txt')
read(1,*)n_ij_n
allocate(ij_n(n_ij_n,2))
read(1,*)
do k=1,n_ij_n
  read(1,*)x,y,ij_n(k,:)
end do
close(1)
! <--- set cv & flux
!
!open(10,file='input/bound.dat')
!do i=0,qnum ! i=0: initial
!  read(10,*)qin0(i),cin0(i)
!end do
!qin0(qnum+1)=qin0(qnum)
!cin0(qnum+1)=cin0(qnum)
!close(10)
!
! ---> set vol from 1d
open(10,file='input/streamConfiguration_inRR.txt')
read(10,*)iend_1d
allocate(ij_1d2d(iend_1d,2,1),xy_1d2d(iend_1d,2))
read(10,*)
do i=1,iend_1d
  read(10,*)itmp,itmp,itmp,itmp,itmp,&
            tmp,tmp,tmp,tmp,tmp,tmp,itmp,&
            xy_1d2d(i,1),xy_1d2d(i,2),ij_1d2d(i,1,1),ij_1d2d(i,2,1)
end do
close(10)
deallocate(ij_1d2d)
!
open(10,file='input/streamFloodplainConnection.txt')
read(10,*)itr_vin
read(10,*)
itmp=0
io=0
do while(io==0)
  read(10,*,iostat=io)
  itmp=itmp+1
end do
close(10)
itmp=itmp-1
n_ivin=itmp
if(mod(n_ivin,itr_vin)/=0)then
  write(*,*)'n%itr/=0 <<<STOP>>>'
  stop
else
n_ivin=n_ivin/itr_vin
n_ivin0=n_ivin
end if
!
allocate(ij_1d2d(n_ivin,2,itr_vin),i_1d(n_ivin,itr_vin))
open(10,file='input/streamFloodplainConnection.txt')
read(10,*)
read(10,*)
do i=1,n_ivin
  do j=1,itr_vin
    read(10,*)tmp,tmp,ij_1d2d(i,1,j),ij_1d2d(i,2,j),i_1d(i,j)
  end do
end do
close(10)
!
allocate(va_in(iend_1d))
open(10,file='./input/flowVolume_ws_DRtoDF.txt')
do i=1,qnum
  read(10,'(a10,i10,*(E20.10e3))')ctmp,k,va_in(:)
end do
close(10)
!
itmp=0
do i=1,n_ivin
  if(0.d0 < va_in(i_1d(i,1)))then
    itmp=itmp+1
  end if
end do
n_ivin=itmp
!
deallocate(ij_1d2d,i_1d)
allocate(ivin(n_ivin))
allocate(ij_1d2d(n_ivin,2,itr_vin),i_1d(n_ivin,itr_vin))
!
open(10,file='input/streamFloodplainConnection.txt')
read(10,*)
read(10,*)
k=0
do i=1,n_ivin0
  j=1
  read(10,*)tmp,tmp,itmp,itmp,itmp1,itmp2
  backspace(10)
  if(0.d0 < va_in(itmp1))then
    k=k+1
    do j=1,itr_vin
      read(10,*)tmp,tmp,ij_1d2d(k,1,j),ij_1d2d(k,2,j),i_1d(k,j)
    end do
    if(i_1d(k,1)==i_1d(k,itr_vin))then
      ivin(k)=i_1d(k,1)
    else
      write(*,*)'n_vin ? <<<STOP>>>'
      stop
    end if
  else
    do j=1,itr_vin
      read(10,*)
    end do
  end if
end do
deallocate(va_in)
!
allocate(va_in0(0:qnum,iend_1d),vcc_in0(0:qnum,iend_1d),vcf_in0(0:qnum,iend_1d), &
         va_in(iend_1d),vcc_in(iend_1d),vcf_in(iend_1d))
va_in0=0.d0
vcc_in0=0.d0
vcf_in0=0.d0
va_in=0.d0
vcc_in=0.d0
vcf_in=0.d0
open(10,file='./input/flowVolume_ws_DRtoDF.txt')
do i=0,qnum-1 ! i=0: initial
  read(10,'(a10,i10,*(E20.10e3))')ctmp,k,va_in0(i,:)
  if(i==0)then
    va_in(:)=va_in0(i,:)
  else
    va_in0(i,:)=va_in0(i,:)-va_in(:)
  end if
  va_in(:)=va_in(:)+va_in0(i,:)
end do
close(10)
va_in0(qnum,:)=va_in0(qnum-1,:)
va_in0(:,:)=va_in0(:,:)/(cellsize*cellsize*qdt)
!
open(10,file='./input/flowVolume_cs_DRtoDF.txt')
do i=0,qnum-1 ! i=0: initial
  read(10,'(a10,i10,*(E20.10e3))')ctmp,k,vcc_in0(i,:)
  if(i==0)then
    vcc_in(:)=vcc_in0(i,:)
  else
    vcc_in0(i,:)=vcc_in0(i,:)-vcc_in(:)
  end if
  vcc_in(:)=vcc_in(:)+vcc_in0(i,:)
end do
close(10)
vcc_in0(qnum,:)=vcc_in0(qnum-1,:)
vcc_in0(:,:)=vcc_in0(:,:)/(cellsize*cellsize*qdt)
!
open(10,file='./input/flowVolume_fs_DRtoDF.txt')
do i=0,qnum-1 ! i=0: initial
  read(10,'(a10,i10,*(E20.10e3))')ctmp,k,vcf_in0(i,:)
  if(i==0)then
    vcf_in(:)=vcf_in0(i,:)
  else
    vcf_in0(i,:)=vcf_in0(i,:)-vcf_in(:)
  end if
  vcf_in(:)=vcf_in(:)+vcf_in0(i,:)
end do
close(10)
vcf_in0(qnum,:)=vcf_in0(qnum-1,:)
vcf_in0(:,:)=vcf_in0(:,:)/(cellsize*cellsize*qdt)
!
open(10,file='va.txt')
write(10,'((i10),(a10),*(i20))')n_ivin*itr_vin,'i_1D',((i_1d(j,k),k=1,itr_vin),j=1,n_ivin)
write(10,'((i10),(a10),*(i20))')n_ivin*itr_vin,'i_2D',((ij_1d2d(j,1,k),k=1,itr_vin),j=1,n_ivin)
write(10,'((i10),(a10),*(i20))')n_ivin*itr_vin,'j_2D',((ij_1d2d(j,2,k),k=1,itr_vin),j=1,n_ivin)
do i=1,qnum
  write(10,'((a10),(i10),*(E20.10e3))')'time(s)=',i*int(qdt),((va_in0(i,i_1d(j,k))/dble(itr_vin),k=1,itr_vin),j=1,n_ivin)
end do
close(10)
!
open(10,file='vcc.txt')
write(10,'((i10),(a10),*(i20))')n_ivin*itr_vin,'i_1D',((i_1d(j,k),k=1,itr_vin),j=1,n_ivin)
write(10,'((i10),(a10),*(i20))')n_ivin*itr_vin,'i_2D',((ij_1d2d(j,1,k),k=1,itr_vin),j=1,n_ivin)
write(10,'((i10),(a10),*(i20))')n_ivin*itr_vin,'j_2D',((ij_1d2d(j,2,k),k=1,itr_vin),j=1,n_ivin)
do i=1,qnum
  write(10,'((a10),(i10),*(E20.10e3))')'time(s)=',i*int(qdt),((vcc_in0(i,i_1d(j,k))/dble(itr_vin),k=1,itr_vin),j=1,n_ivin)
end do
close(10)
!
open(10,file='vcf.txt')
write(10,'((i10),(a10),*(i20))')n_ivin*itr_vin,'i_1D',((i_1d(j,k),k=1,itr_vin),j=1,n_ivin)
write(10,'((i10),(a10),*(i20))')n_ivin*itr_vin,'i_2D',((ij_1d2d(j,1,k),k=1,itr_vin),j=1,n_ivin)
write(10,'((i10),(a10),*(i20))')n_ivin*itr_vin,'j_2D',((ij_1d2d(j,2,k),k=1,itr_vin),j=1,n_ivin)
do i=1,qnum
  write(10,'((a10),(i10),*(E20.10e3))')'time(s)=',i*int(qdt),((vcf_in0(i,i_1d(j,k))/dble(itr_vin),k=1,itr_vin),j=1,n_ivin)
end do
close(10)
!
open(10,file='connection_mod.txt')
write(10,*)n_ivin*itr_vin
do i=1,n_ivin
  do j=1,itr_vin
    x=xllcorner+cellsize*(0.5+ij_1d2d(i,1,j)-1)
    y=yllcorner+cellsize*(0.5+ij_1d2d(i,2,j)-1)
    write(10,'(2(f15.3,'',''),5(i5,'',''))')x,y,(i-1)*itr_vin+j,ij_1d2d(i,1,j),ij_1d2d(i,2,j),i_1d(i,j),ivin(i)
  end do
end do
close(10)
! <--- set vol from 1d
!
! ---> reconfiguration of input vol
deallocate(ij_1d2d)
open(10,file='connection_mod.txt')
read(10,*)n_ivin
allocate(ij_1d2d(n_ivin,2,1))
do i=1,n_ivin
  read(10,*)x,y,itmp,ij_1d2d(i,1,1),ij_1d2d(i,2,1)
end do
!
deallocate(va_in0,vcc_in0,vcf_in0,va_in,vcc_in,vcf_in)
allocate(va_in0(0:qnum,n_ivin),vcc_in0(0:qnum,n_ivin),vcf_in0(0:qnum,n_ivin), &
         va_in(n_ivin),vcc_in(n_ivin),vcf_in(n_ivin))
va_in0=0.d0
vcc_in0=0.d0
vcf_in0=0.d0
va_in=0.d0
vcc_in=0.d0
vcf_in=0.d0
open(10,file='va.txt')
read(10,*)
read(10,*)
read(10,*)
do i=1,qnum ! i=0: initial
  read(10,'((a10),(i10),*(E20.10e3))')ctmp,k,va_in0(i,1:n_ivin)
end do
close(10)
!
open(10,file='vcc.txt')
read(10,*)
read(10,*)
read(10,*)
do i=1,qnum ! i=0: initial
  read(10,'((a10),(i10),*(E20.10e3))')ctmp,k,vcc_in0(i,1:n_ivin)
end do
close(10)
!
open(10,file='vcf.txt')
read(10,*)
read(10,*)
read(10,*)
do i=1,qnum ! i=0: initial
  read(10,'((a10),(i10),*(E20.10e3))')ctmp,k,vcf_in0(i,1:n_ivin)
end do
close(10)
! <--- reconfiguration of input vol
!
! ---> set rain
open(11,file=rname)
read(11,*)n_rpnt
allocate(rseries(iend,jend,0:rnum),rmap(iend,jend,0:n_rpnt),r(iend,jend), &
         x_rpnt(n_rpnt),y_rpnt(n_rpnt),r_rpnt(n_rpnt),d_idw(n_rpnt))
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
fname=trim(adjustl(out_dir))//'rain_map.asc'
cfmt='(*(f10.3))'
call w_fasc_fmt(20,fname,rmap(:,:,0),cfmt,iend,jend,xllcorner,yllcorner,cellsize,-9999.d0)

do k=0,rnum-1
  rseries(:,:,k)=0.0d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  read(11,*)ctmp,r_rpnt(1:n_rpnt)
  do i=1,iend
    do j=1,jend
      if(ibasin(i,j)/=0)then
        rain_tmp=sum(r_rpnt(1:n_rpnt)*rmap(i,j,1:n_rpnt))
        !write(*,'(2i4,10f7.3)')i,j,rain_tmp,r_rpnt(:),rmap(i,j,:)
        rseries(i,j,k)=rain_tmp/rdt/1000.0d0
      end if
    end do
  end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !rseries(1:iend,1:jend,k)=rain_tmp/rdt/1000.0d0
  !where(ibasin(1:iend,1:jend)==0) rseries(1:iend,1:jend,k)=0.0d0
  !write(fname,'(i10.10)')nint(k*rdt)
  !fname='output_rain/rain_'//trim(adjustl(fname))//'.asc'
  !write(*,*)trim(adjustl(fname)),' mm/h'
  !call w_fasc(20,fname,cfmt,rseries(:,:,k)*3600.*1000.0d0,iend,jend,xllcorner,yllcorner,dx,-9999.0d0)
  !open(10,file=fname)
  !do i=1,6
  !  read(10,*)
  !end do
  !do j=jend,1,-1
  !  read(10,*)rseries(1:iend,j,k)
  !  rseries(1:iend,j,k)=rseries(1:iend,j,k)/rdt/1000.0d0
  !end do
end do

fname=trim(adjustl(out_dir))//'rain_sum.asc'
write(*,*)trim(adjustl(fname)),' mm'
r=0.d0
do k=0,rnum-1
 r(:,:)=r(:,:)+(rseries(:,:,k))*rdt*1000.
end do
!write(*,*)'r_ave=',sum(r(:,:))*dx*dx/basin_area
call w_fasc_fmt(20,fname,r(:,:),cfmt,iend,jend,xllcorner,yllcorner,cellsize,-9999.0d0)
r=0.d0
close(11)
rseries(:,:,rnum)=rseries(:,:,rnum-1)
! <--- set rain
!
! ---> initial conditions
h(:,:)=0.d0
u(:,:)=0.d0
v(:,:)=0.d0
uh(:,:)=0.d0
vh(:,:)=0.d0
ho(:,:)=h(:,:)
uo(:,:)=u(:,:)
vo(:,:)=v(:,:)
zo(:,:)=z(:,:)
cc(:,:)=0.d0
cf(:,:)=0.d0
cco(:,:)=cc(:,:)
cfo(:,:)=cf(:,:)
e(:,:)=0.d0
q_a(:,:,1:2)=0.d0
qcc_a(:,:,1:2)=0.d0
qcf_a(:,:,1:2)=0.d0
!hs(:,:)=0.d0
coh(:,:)=coh_tmp
sf(:,:)=0.d0
hmax=0.d0
hsmax=0.d0
umax=0.d0
vmax=0.d0
i_umax=0
j_umax=0
h_umax=0.d0
i_vmax=0
j_vmax=0
h_vmax=0.d0
va_acc=0.d0
vcc_acc=0.d0
vcf_acc=0.d0
q_out=0.d0
qcc_out=0.d0
qcf_out=0.d0
! <--- initial conditions
!
! ---> cal start
dt0=dt
dtmin=dt
time=0.d0
time0=0.d0
do itsp=1,itsp_max
  time=time0+dt
  k=itsp/itsp_q
  l=itsp/itsp_r
  qin=qin0(k)+(qin0(k+1)-qin0(k))*(itsp-k*itsp_q)/dble(itsp_q)
  cin=cin0(k)+(cin0(k+1)-cin0(k))*(itsp-k*itsp_q)/dble(itsp_q)
  r(:,:)=rseries(:,:,l+1)
  va_in(:)=va_in0(k,:)+(va_in0(k+1,:)-va_in0(k,:))*(itsp-k*itsp_q)/dble(itsp_q)
  vcc_in(:)=vcc_in0(k,:)+(vcc_in0(k+1,:)-vcc_in0(k,:))*(itsp-k*itsp_q)/dble(itsp_q)
  vcf_in(:)=vcf_in0(k,:)+(vcf_in0(k+1,:)-vcf_in0(k,:))*(itsp-k*itsp_q)/dble(itsp_q)
  !
  time1=time0
  itr=1
  flxd=0.d0
  !
  call cal_qs(n_ij_u,ij_u,iend,jend,-1, 0,dx,ks,zs,hs,qsx,usmax)
  call cal_qs(n_ij_v,ij_v,iend,jend, 0,-1,dx,ks,zs,hs,qsy,vsmax)
  !flxd=max(usmax*dt,vsmax*dt,ks*dt)
  !flxd=max(flxd, maxval(abs(qsx))*dt, maxval(abs(qsy))*dt)
  !
  call cal_rho(n_ij_cv,ij_cv,iend,jend,sig,rhow,cc,cf,rho,rhom) ! i=i to iend
  !u
  call cal_u(n_ij_u,ij_u,iend,jend,dtlim,dx,g,zo,h,ho,tant_e,hlim,rn,u,uo,uh, &
             rho,sig,dm,csta,tanp,cco,ijofs_uv(1:8,0:4,1),vo,umax,i_umax,j_umax,h_umax)
  !v
  call cal_u(n_ij_v,ij_v,iend,jend,dtlim,dx,g,zo,h,ho,tant_e,hlim,rn,v,vo,vh, &
             rho,sig,dm,csta,tanp,cco,ijofs_uv(1:8,0:4,2),uo,vmax,i_vmax,j_vmax,h_vmax)
  !
  !flxd=max(flxd, umax*h_umax*dt, vmax*h_vmax*dt)
  flxd=max(umax*dt, vmax*dt)
  !tmp=u2*maxval(va_in(:))*(2.d0*g*maxval(va_in(:)))**0.5d0
  tmp=(g*maxval(va_in(:)))**0.5d0
  flxd=max(flxd, tmp*dt)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if ((dx*cfl) < flxd)then
    itr=min(ceiling(flxd/(dx*cfl)),nint(dt0/dtlim))
    if(itr<=0)then
      write(*,*)'itr',itr,flxd,dx*cfl,nint(dt0/dtlim)
      write(*,*)tmp,h_umax,h_vmax,maxval(va_in(:))
      read(*,*)
    end if
    dt=dt0/dble(itr)
    !write(*,'(a,f10.7)', advance='yes')'dt=',dt
    if(dt<dtmin)dtmin=dt
  end if
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  do while (itr>0)
    !
    call cal_hs(n_ij_cv,ij_cv,iend,jend,dt,dx,hs,r(:,:),qsx,qsy,csta)
    where(d(:,:) < hs(:,:))
      h(:,:)=h(:,:)+(hs(:,:)-d(:,:))*(1.d0-csta)
      hs(:,:)=d(:,:)
    !elsewhere
    end where
    !
    call cal_stability(n_ij_cv,ij_cv,iend,jend, &
                       grad_w,h,hs, &
                       (1.d0-csta),rhow,sig,d,tanp,coh,g,sf)
    !where(1.d0<sf(:) .and. d(:)==hs(:))
    !elsewhere
    !  e(:,:)=0.d0
    !end where
    !
    call cal_h(n_ij_cv,ij_cv,iend,jend,dt,dx,h,ho,qx,qy,e,csta,cstad)
    call cal_cc(n_ij_cv,ij_cv,iend,jend,dt,dx,hlim,h,ho,cc,cco,qccx,qccy,e,pc)
    call cal_cf(n_ij_cv,ij_cv,iend,jend,dt,dx,hlim,h,ho,cf,cfo,cc,cco,qcfx,qcfy,e,pf,cstad)
    call cal_z(n_ij_cv,ij_cv,iend,jend,dt,csta,z,zs,e,tant_e)
    !
    ! ---> from 1d
    do k=1,n_ivin
        i=ij_1d2d(k,1,1)
        j=ij_1d2d(k,2,1)
        ho(i,j)=h(i,j)
        cco(i,j)=cc(i,j)
        cfo(i,j)=cf(i,j)
        h(i,j)=ho(i,j)+va_in(k)*dt
        if(0.d0<h(i,j))then
          cc(i,j)=(cco(i,j)       *ho(i,j)          + vcc_in(k)*dt) /                 h(i,j)
          cf(i,j)=((1.d0-cco(i,j))*ho(i,j)*cfo(i,j) + vcf_in(k)*dt) / ((1.d0-cc(i,j))*h(i,j))
        else
          cc(i,j)=0.d0
          cf(i,j)=0.d0
        end if
        va_acc = va_acc+ va_in(k)*dt*cellsize*cellsize
        vcc_acc=vcc_acc+vcc_in(k)*dt*cellsize*cellsize
        vcf_acc=vcf_acc+vcf_in(k)*dt*cellsize*cellsize
    end do
    ! <--- from 1d
    !
    call cal_qs(n_ij_u,ij_u,iend,jend,-1, 0,dx,ks,zs,hs,qsx,usmax)
    call cal_qs(n_ij_v,ij_v,iend,jend, 0,-1,dx,ks,zs,hs,qsy,vsmax)
    !u
    call cal_u(n_ij_u,ij_u,iend,jend,dt,dx,g,zo,h,ho,tant_e,hlim,rn,u,uo,uh, &
               rho,sig,dm,csta,tanp,cco,ijofs_uv(1:8,0:4,1),vo,umax,i_umax,j_umax,h_umax)
    !v
    call cal_u(n_ij_v,ij_v,iend,jend,dt,dx,g,zo,h,ho,tant_e,hlim,rn,v,vo,vh, &
               rho,sig,dm,csta,tanp,cco,ijofs_uv(1:8,0:4,2),uo,vmax,i_vmax,j_vmax,h_vmax)
    !
    ua(1:iend-1,:)=(u(1:iend-1,:)+u(2:iend,:))*0.5d0
    ua(iend,:)=u(iend,:)
    va(:,1:jend-1)=(v(:,1:jend-1)+v(:,2:jend))*0.5d0
    va(:,jend)=v(:,jend)
    uv=(ua*ua+va*va)**0.5d0
    !
    call cal_q(n_ij_u,ij_u,iend,jend,-1, 0,h,u,qx,cc,cf,qccx,qcfx,uh)
    call cal_q(n_ij_v,ij_v,iend,jend, 0,-1,h,v,qy,cc,cf,qccy,qcfy,vh)
    !
    call cal_rho(n_ij_cv,ij_cv,iend,jend,sig,rhow,cc,cf,rho,rhom) ! i=i to iend
    !
    call cal_grad(n_ij_cv,ij_cv,iend,jend,dx,z,grad_z)
    call cal_grad(n_ij_cv,ij_cv,iend,jend,dx,zs+hs+h,grad_w)
    !
    call cal_e(n_ij_cv,ij_cv,iend,jend,dt,dx,z,zs,sig,rho,csta,tanp,cc,cf,e,h,hlim, &
               ua,va,qx,qy,qccx,qccy,qcfx,qcfy,ifld,tant_e,uv) ! output: i=1 to iend-1
    !
    uo(:,:)=u(:,:)
    vo(:,:)=v(:,:)
    ho(:,:)=h(:,:)
    zo(:,:)=z(:,:)
    cco(:,:)=cc(:,:)
    cfo(:,:)=cf(:,:)
    !
    q_a(:,:,1)=q_a(:,:,1)+qx(:,:)*dt
    q_a(:,:,2)=q_a(:,:,2)+qy(:,:)*dt
    qcc_a(:,:,1)=qcc_a(:,:,1)+qccx(:,:)*dt
    qcc_a(:,:,2)=qcc_a(:,:,2)+qccy(:,:)*dt
    qcf_a(:,:,1)=qcf_a(:,:,1)+qcfx(:,:)*dt
    qcf_a(:,:,2)=qcf_a(:,:,2)+qcfy(:,:)*dt
    !
    do k=1,n_ij_w
      q_out=q_out-qx(ij_w(k,1),ij_w(k,2))*dt*dx
      qcc_out=qcc_out-qccx(ij_w(k,1),ij_w(k,2))*dt*dx
      qcf_out=qcf_out-qcfx(ij_w(k,1),ij_w(k,2))*dt*dx
    end do
    do k=1,n_ij_e
      q_out=q_out+qx(ij_e(k,1),ij_e(k,2))*dt*dx
      qcc_out=qcc_out+qccx(ij_e(k,1),ij_e(k,2))*dt*dx
      qcf_out=qcf_out+qcfx(ij_e(k,1),ij_e(k,2))*dt*dx
    end do
    do k=1,n_ij_s
      q_out=q_out-qy(ij_s(k,1),ij_s(k,2))*dt*dx
      qcc_out=qcc_out-qccy(ij_s(k,1),ij_s(k,2))*dt*dx
      qcf_out=qcf_out-qcfy(ij_s(k,1),ij_s(k,2))*dt*dx
    end do
    do k=1,n_ij_n
      q_out  =q_out  +qy(  ij_n(k,1),ij_n(k,2))*dt*dx
      qcc_out=qcc_out+qccy(ij_n(k,1),ij_n(k,2))*dt*dx
      qcf_out=qcf_out+qcfy(ij_n(k,1),ij_n(k,2))*dt*dx
    end do
    !
    where(hmax(:,:)<h(:,:))hmax(:,:)=h(:,:)
    where(hsmax(:,:)<hs(:,:))hsmax(:,:)=hs(:,:)
    !
    time1=time1+dt
    itr=itr-1
    !write(*,'(a,3f10.5,i5)')'in do while1, time,time0,time1',time,time0,time1,itr
  end do ! do while
  !write(*,'(a,3f10.5)')'af do while1, time,time0,time1',time,time0,time1
  !
  time0=time
  dt=dt0
  ! ---> output
  if (mod(itsp, itsp_dout) ==0) then
    ! ---> to disp
    write(*,'(a,f10.3)')      'time(day) =',time/86400.d0
    write(*,'(a,f10.3)')      'time(hour)=',time/3600.d0
    write(*,'(a,i10,2(a,f7.3))')'time(sec) =',nint(time),', rain(mm/h)=',maxval(r(:,:)*1000*3600.),&
    ', va_in(m/s)=',maxval(va_in(:)/dble(itr_vin))
    write(*,'(2(a,f10.3),a,2i5)')'umax=',umax,', h@umax=',h_umax,', i,j=',i_umax,j_umax
    write(*,'(2(a,f10.3),a,2i5)')'vmax=',vmax,', h@vmax=',h_vmax,', i,j=',i_vmax,j_vmax
    write(*,'(2(a,f10.5))')'dt_min=',dtmin,', dt_lim=',dtlim
    write(*,'(a,3e20.5e3)')'va_in, vcc_in, vcf_in',va_acc,vcc_acc,vcf_acc
    write(*,'(a,3e20.5e3)')'vasum, vccsum, vcfsum',q_out+sum(h)*cellsize*cellsize,&
    qcc_out+sum(h*cc+(z-zini)*csta)*cellsize*cellsize,&
    qcf_out+sum(h*(1.d0-cc)*cf)*cellsize*cellsize
    write(*,'(a,3e20.5e3)')'q_out,qcc_out,qcf_out',q_out,qcc_out,qcf_out
    ! <--- to disp
    umax=0.d0
    vmax=0.d0
    dtmin=dt
    !
    ! ---> flux at section
    cline=''
    ! ---> sec i
    do i=1,nsci
      j=(i-1)*3
      mvals(j+1)=sum(q_a(  iffx(i,1),iffx(i,2):iffx(i,3),1))*dx
      mvals(j+2)=sum(qcc_a(iffx(i,1),iffx(i,2):iffx(i,3),1))*dx
      mvals(j+3)=sum(qcf_a(iffx(i,1),iffx(i,2):iffx(i,3),1))*dx
    end do
    ! <--- sec i
    ! ---> sec j
    do i=1,nscj
      j=nsci*3+(i-1)*3
      mvals(j+1)=sum(q_a(  iffy(i,2):iffy(i,3),iffy(i,1),2))*dx
      mvals(j+2)=sum(qcc_a(iffy(i,2):iffy(i,3),iffy(i,1),2))*dx
      mvals(j+3)=sum(qcf_a(iffy(i,2):iffy(i,3),iffy(i,1),2))*dx
    end do
    ! <--- sec i
    write(30,'(a,i5,20f15.3)')'time=',nint(time),mvals
    ! <--- flux at section
    !
    ! ---> to files
    if(mod(itsp, itsp_fout) ==0) then
      ua(1:iend-1,:)=(u(1:iend-1,:)+u(2:iend,:))*0.5d0
      ua(iend,:)=u(iend,:)
      va(:,1:jend-1)=(v(:,1:jend-1)+v(:,2:jend))*0.5d0
      va(:,jend)=v(:,jend)
      uv=(ua*ua+va*va)**0.5d0
      cfmt='(*(f10.3))'
      write(ctmp,'(i10.10)')idnint(time)
      fname=trim(adjustl(out_dir))//'res_h_'//trim(adjustl(ctmp))//'.asc'
      call w_fasc_fmt(20,fname,h(:,:),cfmt,iend,jend,xllcorner,yllcorner,cellsize,-9999.d0)
      fname=trim(adjustl(out_dir))//'res_hs_'//trim(adjustl(ctmp))//'.asc'
      !call w_fasc_fmt(20,fname,hs(:,:),cfmt,iend,jend,xllcorner,yllcorner,cellsize,-9999.d0)
      fname=trim(adjustl(out_dir))//'res_hmax.asc'
      call w_fasc_fmt(20,fname,hmax(:,:),cfmt,iend,jend,xllcorner,yllcorner,cellsize,-9999.d0)
      fname=trim(adjustl(out_dir))//'res_hsmax.asc'
      call w_fasc_fmt(20,fname,hsmax(:,:),cfmt,iend,jend,xllcorner,yllcorner,cellsize,-9999.d0)
      fname=trim(adjustl(out_dir))//'res_sf_'//trim(adjustl(ctmp))//'.asc'
      !call w_fasc_fmt(20,fname,sf(:,:),cfmt,iend,jend,xllcorner,yllcorner,cellsize,-9999.d0)
      fname=trim(adjustl(out_dir))//'res_cc_'//trim(adjustl(ctmp))//'.asc'
      call w_fasc_fmt(20,fname,cc(:,:),cfmt,iend,jend,xllcorner,yllcorner,cellsize,-9999.d0)
      fname=trim(adjustl(out_dir))//'res_cf_'//trim(adjustl(ctmp))//'.asc'
      !call w_fasc_fmt(20,fname,cf(:,:),cfmt,iend,jend,xllcorner,yllcorner,cellsize,-9999.d0)
      fname=trim(adjustl(out_dir))//'res_dz_'//trim(adjustl(ctmp))//'.asc'
      call w_fasc_fmt(20,fname,z(:,:)-zini(:,:),cfmt,iend,jend,xllcorner,yllcorner,cellsize,-9999.d0)
      !fname=trim(adjustl(out_dir))//'res_e_'//trim(adjustl(ctmp))//'.asc'
      !call w_fasc_fmt(20,fname,e(:,:),'(*(f10.3))',iend,jend,xllcorner,yllcorner,cellsize,-9999.d0)
      fname=trim(adjustl(out_dir))//'res_qx_a_'//trim(adjustl(ctmp))//'.asc'
      !call w_fasc_fmt(20,fname,q_a(:,:,1),cfmt,iend,jend,xllcorner-dx*0.5d0,yllcorner,cellsize,-9999.d0)
      fname=trim(adjustl(out_dir))//'res_qy_a_'//trim(adjustl(ctmp))//'.asc'
      !call w_fasc_fmt(20,fname,q_a(:,:,2),cfmt,iend,jend,xllcorner,yllcorner-dx*0.5d0,cellsize,-9999.d0)
      fname=trim(adjustl(out_dir))//'res_qccx_a_'//trim(adjustl(ctmp))//'.asc'
      !call w_fasc_fmt(20,fname,qcc_a(:,:,1),cfmt,iend,jend,xllcorner-dx*0.5d0,yllcorner,cellsize,-9999.d0)
      fname=trim(adjustl(out_dir))//'res_qccy_a_'//trim(adjustl(ctmp))//'.asc'
      !call w_fasc_fmt(20,fname,qcc_a(:,:,2),cfmt,iend,jend,xllcorner,yllcorner-dx*0.5d0,cellsize,-9999.d0)
      fname=trim(adjustl(out_dir))//'res_qcfx_a_'//trim(adjustl(ctmp))//'.asc'
      !call w_fasc_fmt(20,fname,qcf_a(:,:,1),cfmt,iend,jend,xllcorner-dx*0.5d0,yllcorner,cellsize,-9999.d0)
      fname=trim(adjustl(out_dir))//'res_qcfy_a_'//trim(adjustl(ctmp))//'.asc'
      !call w_fasc_fmt(20,fname,qcf_a(:,:,2),cfmt,iend,jend,xllcorner,yllcorner-dx*0.5d0,cellsize,-9999.d0)
      !fname=trim(adjustl(out_dir))//'res_uv_'//trim(adjustl(ctmp))//'.asc'
      !call w_fasc_fmt(20,fname,uv(:,:),cfmt,iend,jend,xllcorner,yllcorner,cellsize,-9999.d0)
      ! q_a(:,:,1:2)=0.d0
      ! qcc_a(:,:,1:2)=0.d0
      ! qcf_a(:,:,1:2)=0.d0
    end if
    ! <--- to files
  end if
  ! <--- output
end do
! <--- cal start
close(30)
write(*,*)'--- nomal end ---'
!
stop
end
!
!
!
subroutine cal_u(n_1d,ij_1d,iend,jend,dt,dx,g,zo,h,ho,tan_a,hlim,rn,u,uo,uh, &
                 rho,sig,dm,csta,tanp,cco,ijofs,v,umax,i_umax,j_umax,h_umax)
! FX q1   q2   q3    qe
!    ->   ->   ->    ->
! CV | z1 | z2 | ....| ze |
use omp_lib
implicit none
real(8), parameter :: f1d5=1.d0/5.d0,f1d3=1.d0/3.d0,f5d3=5.d0/3.d0,f2d3=2.d0/3.d0, &
                      kf=0.16d0,kd=0.0828d0,e2=0.85d0**2.d0
integer :: ni,n_1d,ij_1d(n_1d,2),i,iend,j,jend,i_umax,j_umax,ijofs(8,0:4),io(4),jo(4),i2,j2,iloc(2),icase
real(8) :: va,dt,dx,g,x1,x2,x3,hlim,rn,umax,h_umax,htmp(iend,jend), &
           zo(iend,jend),h(iend,jend),ho(iend,jend),tan_a(iend,jend), &
           u(iend,jend),uo(iend,jend),v(iend,jend),h1,hnw,c1,rho1,z1,z2,tan1, &
           cco(iend,jend),rho(iend,jend),uh(iend,jend), &!,vh(iend,jend)
           sig,csta,dm,taub,tauy,fb,fb2,tanp,fd,ff,cost,rhom,dhdx
!
umax=0.d0
u(:,:)=0.d0
htmp(:,:)=0.d0
h_umax=0.d0
i_umax=0
j_umax=0
!$omp parallel do private(ni,i,j,i2,j2,icase,io,jo, &
!$omp                     va,dhdx,h1,hnw,c1,rho1,z1,z2,tan1,x1,x2,x3, &
!$omp                     cost,tauy,fd,ff,fb2,fb,taub,rhom)
do ni=1,n_1d
  i=ij_1d(ni,1)
  j=ij_1d(ni,2)
  i2=i+ijofs(4,0)
  j2=j+ijofs(8,0)
  if(ho(i,j)<hlim .and. ho(i2,j2)<hlim)then ! .or. &
    !h(i,j)<hlim .and. h(i2,j2)<hlim)then
    u(i,j)=0.d0
    cycle
  end if
  if(0.d0<uo(i,j) .or. uo(i,j)==0.d0 .and. zo(i,j)+ho(i,j) <= zo(i2,j2) +ho(i2,j2) )then
    h1=ho(i2,j2)
    hnw=h(i2,j2)
    c1=cco(i2,j2)
    rho1=rho(i2,j2)
    z1=zo(i2,j2)
    z2=zo(i,j)
    tan1=tan_a(i,j)
  else if(uo(i,j)<0.d0 .or. uo(i,j)==0.d0 .and. zo(i2,j2) +ho(i2,j2) < zo(i,j)+ho(i,j))then
    h1=ho(i,j)
    hnw=h(i,j)
    c1=cco(i,j)
    rho1=rho(i,j)
    z1=zo(i,j)
    z2=zo(i2,j2)
    tan1=tan_a(i2,j2)
  else
    write(*,*)'@cal_u',uo(i,j),h1,i,j
    read(*,*)
  end if
  if(h1<hlim .or. z1+h1<=z2)then! .or. hnw<hlim .or. z1+hnw<=z2)then
    u(i,j)=0.d0
    cycle
  end if
  h1=(ho(i,j)+ho(i2,j2))/2.d0
  hnw=(h(i,j)+h(i2,j2))/2.d0
  c1=(cco(i,j)+cco(i2,j2))/2.d0
  rho1=(rho(i,j)+rho(i2,j2))/2.d0
  tan1=(tan_a(i,j)+tan_a(i2,j2))/2.d0
  !if(h1<hlim .or. z1+h1<=z2)then! .or. hnw<hlim .or. z1+hnw<=z2)then
  !  u(i,j)=0.d0
  !  cycle
  !end if
  io(1:4)=ijofs(1:4,0)
  jo(1:4)=ijofs(5:8,0)
  va=v(i+io(1),j+jo(1))+ &
     v(i+io(2),j+jo(2))+ &
     v(i+io(3),j+jo(3))+ &
     v(i+io(4),j+jo(4)) / 4.d0
  if(0<uo(i,j))then
    if(0<va)then !case1
      icase=1
      io(1:4)=ijofs(1:4,icase)
      jo(1:4)=ijofs(5:8,icase)
    else         !case2
      icase=2
      io(1:4)=ijofs(1:4,icase)
      jo(1:4)=ijofs(5:8,icase)
    end if
  else
    if(0<va)then !case3
      icase=3
      io(1:4)=ijofs(1:4,icase)
      jo(1:4)=ijofs(5:8,icase)
    else         !case4
      icase=4
      io(1:4)=ijofs(1:4,icase)
      jo(1:4)=ijofs(5:8,icase)
    end if
  end if
  !@u (i1,j1),(i2,j2),(i3,j3),(i4,j4) | @v (i1,j1),(i2,j2),(i3,j3),(i4,j4)
  !va:( 0, 0),( 0, 1),(-1, 1),(-1, 0) | ua:( 0, 0),( 1, 0),( 1,-1),( 0,-1)
  ! 1:( 0, 0),(-1, 0),( 0, 0),( 0,-1) |  1:( 0, 0),( 0,-1),( 0, 0),(-1, 0)
  ! 2:( 0, 0),(-1, 0),( 0, 1),( 0, 0) |  2:( 0, 0),( 0,-1),( 1, 0),( 0, 0)
  ! 3:( 1, 0),( 0, 0),( 0, 0),( 0,-1) |  3:( 0, 1),( 0, 0),( 0, 0),(-1, 0)
  ! 4:( 1, 0),( 0, 0),( 0, 1),( 0, 0) |  4:( 0, 1),( 0, 0),( 1, 0),( 0, 0)
  !for u
  !va=v(i,j)+v(i,j+1)+v(i-1,j+1)+v(i-1,j)
  !0<u,0<va: u(i,j)*(m(i  ,j)-m(i-1,j))/dx + va(i,j)*(m(i,  j)-m(i,j-1))/dy
  !0<u,va<0: u(i,j)*(m(i  ,j)-m(i-1,j))/dx + va(i,j)*(m(i,j+1)-m(i,  j))/dy
  !u<0,0<va: u(i,j)*(m(i+1,j)-m(i  ,j))/dx + va(i,j)*(m(i,  j)-m(i,j-1))/dy
  !u<0,va<0: u(i,j)*(m(i+1,j)-m(i  ,j))/dx + va(i,j)*(m(i,j+1)-m(i,  j))/dy
  !for v
  !ua=u(i,j)+u(i+1,j)+u(i+1,j-1)+u(i,j-1)
  !0<v,0<ua: v(i,j)*(n(i,  j)-n(i,j-1))/dy + ua(i,j)*(n(i  ,j)-n(i-1,j))/dx
  !0<v,ua<0: v(i,j)*(n(i,  j)-n(i,j-1))/dy + ua(i,j)*(n(i+1,j)-n(i  ,j))/dx
  !v<0,0<ua: v(i,j)*(n(i,j+1)-n(i,  j))/dy + ua(i,j)*(n(i  ,j)-n(i-1,j))/dx
  !v<0,ua<0: v(i,j)*(n(i,j+1)-n(i,  j))/dy + ua(i,j)*(n(i+1,j)-n(i  ,j))/dx
  x1=uo(i,j)*(uh(i+io(1),j+jo(1))-uh(i+io(2),j+jo(2)))/dx + &
          va*(uh(i+io(3),j+jo(3))-uh(i+io(4),j+jo(4)))/dx
  !
  dhdx=(zo(i,j)+ho(i,j)-zo(i2,j2)-ho(i2,j2))/dx
  x2=g*h1*dhdx
  !
  cost=cos(atan(tan1))
  tauy=(c1/csta)**f1d5*(sig-rho1)*c1*g*h1*cost*tanp
  !if(c1<0.05d0)then
  !  x3=-g*rn*rn*uo(i,j)*abs(uo(i,j))/h1**f1d3
  !else
  if(c1<0.05d0)then
    fb=(6.d0+2.5d0*log(h1/dm))**(-2.d0)
  else
    fd=kd*sig/rho1*(1.d0-e2)*c1**(f1d3)
    ff=kf*(1.d0-c1)**f5d3/c1**(f2d3)
    fb=6.25d0*(fd+ff)*(h1/dm)**(-2.d0)
    !fb=6.25d0*(fd+ff)*(max(h1,dm)/dm)**(-2.d0)
    fb2=(6.d0+2.5d0*log(h1/dm))**(-2.d0)
    if(fb<fb2)fb=fb2
  end if
  taub=tauy+rho1*fb*(uo(i,j)*uo(i,j)+va*va)
  rhom=(sig-rho1)*c1+rho1
  !if(uo(i,j)*uo(i,j)+va*va<0.00001d0)then
  if(uo(i,j)*uo(i,j)+va*va<g*hlim)then
    x3=0.d0
  else
    !x3=taub*uo(i,j)/(uo(i,j)*uo(i,j)+va*va)**0.5d0/rhom
    x3=taub/(2.d0*rhom*(uo(i,j)*uo(i,j)+va*va)**0.5d0)
  end if
  !u(i,j)=(uo(i,j)*h1+(-x1-x2-x3)*dt)/hnw
  u(i,j)=(uo(i,j)*(h1-x3*dt)-(x1+x2)*dt)/(hnw+x3*dt)
  !
  if(u(i,j)*uo(i,j)<0.d0)then
    u(i,j)=0.d0
  end if
  htmp(i,j)=hnw
end do
!$omp end parallel do
umax=maxval(abs(u(:,:)))
iloc=maxloc(abs(u(:,:)))
i_umax=iloc(1)
j_umax=iloc(2)
h_umax=htmp(i_umax,j_umax)
!
end subroutine cal_u
!
!
!
subroutine cal_q(n_1d,ij_1d,iend,jend,iofs,jofs,h,u,qx,cc,cf,qccx,qcfx,uh)
use omp_lib
implicit none
integer :: ni,n_1d,ij_1d(n_1d,2),i,iend,j,jend,iofs,jofs
real(8) :: h(iend,jend),u(iend,jend),qx(iend,jend),cc(iend,jend),cf(iend,jend), &
           qccx(iend,jend),qcfx(iend,jend),uh(iend,jend)
!
!u:iofs=-1,jofs=-0
!v:iofs= 0,jofs=-1
!$omp parallel do private(ni,i,j)
do ni=1,n_1d
  i=ij_1d(ni,1)
  j=ij_1d(ni,2)
  uh(i,j)=u(i,j)*(h(i,j)+h(i+iofs,j+jofs))*0.5d0
  if(u(i,j)>0.d0)then
    qx(i,j)=h(i+iofs,j+jofs)*u(i,j)
    qccx(i,j)=qx(i,j)*cc(i+iofs,j+jofs)
    qcfx(i,j)=qx(i,j)*cf(i+iofs,j+jofs)*(1.d0-cc(i+iofs,j+jofs))
  else
    qx(i,j)=h(i,j)*u(i,j)
    qccx(i,j)=qx(i,j)*cc(i,j)
    qcfx(i,j)=qx(i,j)*cf(i,j)*(1.d0-cc(i,j))
  end if
end do
!$omp end parallel do
!
end subroutine cal_q
!
!
!
subroutine cal_qs(n_1d,ij_1d,iend,jend,iofs,jofs,dx,ks,zs,hs,qsx,umax)
use omp_lib
implicit none
integer :: ni,n_1d,ij_1d(n_1d,2),i,iend,j,jend,iofs,jofs,i1,j1
real(8) :: zs(iend,jend),hs(iend,jend),qsx(iend,jend), &
           dx,ks,grad,u,umax,hs1
!
!u:iofs=-1,jofs=-0
!v:iofs= 0,jofs=-1
umax=0.d0
!$omp parallel do private(ni,i,j,i1,j1,grad,hs1,u) reduction(max: umax)
do ni=1,n_1d
  i=ij_1d(ni,1)
  j=ij_1d(ni,2)
  i1=i+iofs
  j1=j+jofs
  grad=(zs(i1,j1)+hs(i1,j1)-(zs(i,j)+hs(i,j)))/dx
  if(grad>0.d0)then
    hs1=hs(i1,j1)
  else
    hs1=hs(i,j)
  end if
  if(hs1<=0.d0)then
    qsx(i,j)=0.d0
    u=0.d0
  else
    qsx(i,j)=hs1*grad*ks
    u=qsx(i,j)/hs1
  end if
  if(umax<abs(u))then
    umax=abs(u)
  end if
end do
!$omp end parallel do
!
end subroutine cal_qs
!
!
!
subroutine cal_hs(n_ij_cv,ij_cv,iend,jend,dt,dx,hs,r,qsx,qsy,csta)
use omp_lib
implicit none
integer :: ni,n_ij_cv,ij_cv(n_ij_cv,2),i,iend,j,jend
real(8) :: dt,dx,csta,hs(iend,jend),qsx(iend,jend),qsy(iend,jend),r(iend,jend)
!
!$omp parallel do private(ni,i,j)
do ni=1,n_ij_cv
  i=ij_cv(ni,1)
  j=ij_cv(ni,2)
  hs(i,j)=hs(i,j) &
          - (qsx(i+1,j)-qsx(i,j))/dx*dt/(1.d0-csta) &
          - (qsy(i,j+1)-qsy(i,j))/dx*dt/(1.d0-csta) &
          + r(i,j)*dt
  if(hs(i,j)<0.d0)then
    write(*,'(a,2i5,E20.10e3)')'hs<0 (i,j,h)',i,j,hs(i,j)
    hs(i,j)=0.d0
  end if
end do
!$omp end parallel do
!
end subroutine cal_hs
!
!
!
subroutine cal_stability(n_ij_cv,ij_cv,iend,jend,grads,h,hs,lamda,rho,sig,d,tanp,c,g,sf)
use omp_lib
implicit none
real(8), parameter :: PI=acos(-1.0d0),deg2rad=PI/180.0d0,rad2deg=180.0d0/PI
integer :: ni,n_ij_cv,ij_cv(n_ij_cv,2),i,j,iend,jend
real(8) :: lamda,rho,sig,tanp,g,grads(iend,jend), &
           h(iend,jend),hs(iend,jend),c(1:iend,1:jend),sf(iend,1:jend),pw(iend,jend), &
           d(iend,jend),grav,resi,cost,csta,sint,tant
!
csta=1.d0-lamda
!$omp parallel do private(ni,i,j,sint,cost,tant,grav,resi)
do ni=1,n_ij_cv
  i=ij_cv(ni,1)
  j=ij_cv(ni,2)
  sint=sin(atan(grads(i,j)))
  cost=cos(atan(grads(i,j)))
  tant=grads(i,j)
  pw(i,j)=0.d0
  !G=rho*g*D*sint*(sig/rho*csta+(1.-hs/D)*pw+hs/D*(1.-csta)+h/D)
  !R=rho*g*D*cost*(sig/rho*csta+(1.-hs/D)*pw-hs/D*csta)*tanp+c
  if(0.d0<d(i,j))then
    grav=rho*g*d(i,j)*sint*(sig/rho*(1.d0-lamda)+(1.d0-hs(i,j)/d(i,j))*pw(i,j)+hs(i,j)/d(i,j)*lamda+h(i,j)/d(i,j))
    resi=rho*g*d(i,j)*cost*(sig/rho*(1.d0-lamda)+(1.d0-hs(i,j)/d(i,j))*pw(i,j)-hs(i,j)/d(i,j)*(1.d0-lamda))*tanp+c(i,j)
    sf(i,j)=grav/resi !G/R
  else
    sf(i,j)=0.d0
  end if
  !
end do
!$omp end parallel do
!
end subroutine cal_stability
!
!
!
subroutine cal_h(n_ij_cv,ij_cv,iend,jend,dt,dx,h,ho,qx,qy,e,csta,cstad)
use omp_lib
implicit none
integer :: ni,n_ij_cv,ij_cv(n_ij_cv,2),i,iend,j,jend
real(8) :: dt,dx,csta,cstad,h(iend,jend),ho(iend,jend),qx(iend,jend),qy(iend,jend),e(iend,jend)
!
!$omp parallel do private(ni,i,j)
do ni=1,n_ij_cv
  i=ij_cv(ni,1)
  j=ij_cv(ni,2)
  if(e(i,j)>=0.d0)then
    h(i,j)=ho(i,j)-(qx(i+1,j)-qx(i,j))/dx*dt-(qy(i,j+1)-qy(i,j))/dx*dt+e(i,j)/csta*dt
  else
    h(i,j)=ho(i,j)-(qx(i+1,j)-qx(i,j))/dx*dt-(qy(i,j+1)-qy(i,j))/dx*dt+e(i,j)/cstad*dt
  end if
  if(h(i,j)<0.d0)then
    write(*,'(a,2i5,4E15.5e3)')'h<0 (i,j,h,e)',i,j,h(i,j),e(i,j),(qx(i+1,j)-qx(i,j)),(qy(i,j+1)-qy(i,j))
    h(i,j)=0.d0
  end if
end do
!$omp end parallel do
!
end subroutine cal_h
!
!
!
subroutine cal_cc(n_ij_cv,ij_cv,iend,jend,dt,dx,hlim,h,ho,cc,cco,qccx,qccy,e,pc)
use omp_lib
implicit none
integer :: ni,n_ij_cv,ij_cv(n_ij_cv,2),i,iend,j,jend
real(8) :: dt,dx,pc,hlim,h(iend,jend),ho(iend,jend),cc(iend,jend),cco(iend,jend), &
           qccx(iend,jend),qccy(iend,jend),e(iend,jend)
!
!$omp parallel do private(ni,i,j)
do ni=1,n_ij_cv
  i=ij_cv(ni,1)
  j=ij_cv(ni,2)
  if(h(i,j)<hlim)then
    cc(i,j)=0.d0
  else
    if(e(i,j)>=0.d0)then
      cc(i,j)=(cco(i,j)*ho(i,j)-(qccx(i+1,j)-qccx(i,j))/dx*dt-(qccy(i,j+1)-qccy(i,j))/dx*dt+pc*e(i,j)*dt)/h(i,j)
    else
      cc(i,j)=(cco(i,j)*ho(i,j)-(qccx(i+1,j)-qccx(i,j))/dx*dt-(qccy(i,j+1)-qccy(i,j))/dx*dt+e(i,j)*dt)/h(i,j)
    end if
    if(cc(i,j)<0.d0)then
      if(cc(i,j)<-0.001d0)write(*,'(a,2i5,4f15.10)')'cc<0 (i,j,cco*ho,cc*h,e,dqx+dqy)', &
                          i,j,cco(i,j)*ho(i,j),cc(i,j)*h(i,j),e(i,j), &
                          (qccx(i+1,j)-qccx(i,j))/dx*dt + (qccy(i,j+1)-qccy(i,j))/dx*dt
      !read(*,*)
      cc(i,j)=0.d0
    end if
  end if
end do
!$omp end parallel do
!
end subroutine cal_cc
!
!
!
subroutine cal_cf(n_ij_cv,ij_cv,iend,jend,dt,dx,hlim,h,ho,cf,cfo,cc,cco,qcfx,qcfy,e,pf,cstad)
use omp_lib
implicit none
integer :: ni,n_ij_cv,ij_cv(n_ij_cv,2),i,iend,j,jend
real(8) :: dt,dx,pf,cstad,hlim,h(iend,jend),ho(iend,jend),cf(iend,jend),cfo(iend,jend), &
           qcfx(iend,jend),qcfy(iend,jend),e(iend,jend),cc(iend,jend),cco(iend,jend)
!
!$omp parallel do private(ni,i,j)
do ni=1,n_ij_cv
  i=ij_cv(ni,1)
  j=ij_cv(ni,2)
  if(h(i,j)<hlim)then
    cf(i,j)=0.d0
  else
    if(e(i,j)>=0.d0)then
      cf(i,j)=(cfo(i,j)*(1.d0-cco(i,j))*ho(i,j)-(qcfx(i+1,j)-qcfx(i,j))/dx*dt-(qcfy(i,j+1)-qcfy(i,j))/dx*dt+pf*e(i,j)*dt) &
            /((1.d0-cc(i,j))*h(i,j))
    else
      cf(i,j)=(cfo(i,j)*(1.d0-cco(i,j))*ho(i,j)-(qcfx(i+1,j)-qcfx(i,j))/dx*dt &
                                               -(qcfy(i,j+1)-qcfy(i,j))/dx*dt+((1.d0-cstad)/cstad)*cfo(i,j)*e(i,j)*dt) &
            /((1.d0-cc(i,j))*h(i,j))
    end if
  end if
  if(cf(i,j)<0.d0)then
    write(*,'(a,2i5,4E20.10e3)')'cc<0 (i,j,cc,cf,h)',i,j,cc(i,j),cf(i,j),e(i,j),h(i,j)
    write(*,'(4f15.10)')qcfx(i+1,j),qcfx(i,j),qcfy(i,j+1),qcfy(i,j)
    !read(*,*)
    cf(i,j)=0.d0
  end if
end do
!$omp end parallel do
!
end subroutine cal_cf
!
!
!
subroutine cal_e(n_ij_cv,ij_cv,iend,jend,dt,dx,z,zs,sig,rho,csta,tanp,cc,cf,e,h,hlim,ua,va, &
                 qx,qy,qccx,qccy,qcfx,qcfy,ifld,tant_e,uv)
use omp_lib
implicit none
real(8), parameter :: PI=acos(-1.d0),deg2rad=PI/180.d0,rad2deg=180.d0/PI
integer :: ni,n_ij_cv,ij_cv(n_ij_cv,2),i,iend,j,jend,ifld(iend,jend)
real(8) :: dt,dx,sig,tanp,tant,tane,csta,emax,dzdx,dzdy,hlim, &
           z(iend,jend),zs(iend,jend),cc(iend,jend),cf(iend,jend),e(iend,jend), &
           h(iend,jend),rho(iend,jend),wl(iend,jend), &
           qx(iend,jend),qy(iend,jend),qccx(iend,jend),qccy(iend,jend),qcfx(iend,jend),qcfy(iend,jend), &
           emin1,emin2,emin3,ua(iend,jend),va(iend,jend),uv(iend,jend),tant_e(iend,jend)
!
wl=z+h
!$omp parallel do private(ni,i,j,dzdx,dzdy,tant,tane, &
!$omp                     emax,emin1,emin2,emin3)
do ni=1,n_ij_cv
  i=ij_cv(ni,1)
  j=ij_cv(ni,2)
  if(h(i,j)<hlim .and. ifld(i,j)==0)then
    e(i,j)=0.d0
    cycle
  end if
  if(ua(i,j)>=0.d0)then
    dzdx=(wl(i+1,j)-wl(i,j))/dx
  else
    dzdx=(wl(i,j)-wl(i-1,j))/dx
  end if
  if(va(i,j)>=0.d0)then
    dzdy=(wl(i,j+1)-wl(i,j))/dx
  else
    dzdy=(wl(i,j)-wl(i,j-1))/dx
  end if
  if(uv(i,j)==0.d0)then
    tant=0.d0
  else
    !tant=max(0.d0,-(dzdx*ua(i,j)+dzdy*va(i,j))/((ua(i,j)*ua(i,j)+va(i,j)*va(i,j))**0.5d0))
    tant=(dzdx*dzdx+dzdy*dzdy)**0.5d0
  end if
  !tant=(dzdx*dzdx+dzdy*dzdy)**0.5d0
  tant_e(i,j)=tant
  tane=(sig/rho(i,j)-1.d0)*cc(i,j)*tanp/((sig/rho(i,j)-1.d0)*cc(i,j)+1.d0)
  !tan(a-b)=(tan(a)-tan(b))/(1+tan(a)*tan(b))
  e(i,j)=csta*(tant-tane)/(1.d0+tant*tane)*uv(i,j)
  emax=(z(i,j)-zs(i,j))/dt*csta*cos(atan(tant))
  if(e(i,j)>=0.d0)then ! ero
    e(i,j)=min(emax,e(i,j))
    !if(ifld(i,j)==1)e(i,j)=emax
  else                 ! depo
    emin1=min(0.d0,(-h(i,j)/dt*csta+(qx(i+1,j)-qx(i,j))/dx+(qy(i,j+1)-qy(i,j))/dx)*csta) ! cstad -> csta
    emin2=min(0.d0,-cc(i,j)*h(i,j)/dt+(qccx(i+1,j)-qccx(i,j))/dx+(qccy(i,j+1)-qccy(i,j))/dx)
    if(cf(i,j)>0.d0)then
      emin3= &
      (-cf(i,j)*(1.d0-cc(i,j))*h(i,j)/dt+(qcfx(i+1,j)-qcfx(i,j))/dx+(qcfy(i,j+1)-qcfy(i,j))/dx) &
      *csta/((1.d0-csta)*cf(i,j))
    else
      emin3=e(i,j)
    end if
    emin3=min(0.d0,emin3)
    e(i,j)=max(emin1,emin2,emin3,e(i,j))
  end if
  if(e(i,j)>emax)then
    write(*,'(a,5f10.3)')'e <--- max',e(i,j),emax,emin1,emin2,emin3
    read(*,*)
    e(i,j)=emax
  end if
end do
!$omp end parallel do
!
end subroutine cal_e
!
!
!
subroutine cal_grad(n_ij_cv,ij_cv,iend,jend,dx,z,grad)
use omp_lib
implicit none
real(8), parameter :: PI=acos(-1.d0),deg2rad=PI/180.d0,rad2deg=180.d0/PI
integer :: ni,n_ij_cv,ij_cv(n_ij_cv,2), i,iend,j,jend
real(8) :: dx,dx2,dzdx,dzdy,z(iend,jend),grad(iend,jend)
!
dx2=dx*2.d0
!do j=2,jend-1
!do i=2,iend-1
!$omp parallel do private(ni,i,j,dzdx,dzdy)
do ni=1,n_ij_cv
  i=ij_cv(ni,1)
  j=ij_cv(ni,2)
  dzdx=(z(i+1,j)-z(i-1,j))/dx2
  dzdy=(z(i,j+1)-z(i,j-1))/dx2
  grad(i,j)=(dzdx*dzdx+dzdy*dzdy)**0.5d0
end do
!$omp end parallel do
!
end subroutine cal_grad
!
!
!
subroutine cal_z(n_ij_cv,ij_cv,iend,jend,dt,csta,z,zs,e,tant_e)
use omp_lib
implicit none
real(8), parameter :: PI=acos(-1.d0),deg2rad=PI/180.d0,rad2deg=180.d0/PI
integer :: ni,n_ij_cv,ij_cv(n_ij_cv,2),i,iend,j,jend
real(8) :: dt,csta,tant,cost,z(iend,jend),zs(iend,jend),e(iend,jend),tant_e(iend,jend)
!
!$omp parallel do private(ni,i,j,tant,cost)
do ni=1,n_ij_cv
  i=ij_cv(ni,1)
  j=ij_cv(ni,2)
  tant=tant_e(i,j)
  cost=cos(atan(tant))
  z(i,j)=z(i,j)-e(i,j)/csta/cost*dt
  if(z(i,j)<zs(i,j))then
    if(z(i,j)-zs(i,j)<-0.0001d0)write(*,'(a,2i5,5f15.5)')'z < zs @ (i,j)=',i,j,e(i,j),z(i,j)-zs(i,j)
    z(i,j)=zs(i,j)
  end if
end do
!$omp end parallel do
!
end subroutine cal_z
!
!
!
subroutine cal_rho(n_ij_cv,ij_cv,iend,jend,sig,rhow,cc,cf,rho,rhom)
use omp_lib
implicit none
integer :: ni,n_ij_cv,ij_cv(n_ij_cv,2),i,iend,j,jend
real(8) :: sig,rhow,cc(iend,jend),cf(iend,jend),rho(iend,jend),rhom(iend,jend)
!
!$omp parallel do private(ni,i,j)
do ni=1,n_ij_cv
  i=ij_cv(ni,1)
  j=ij_cv(ni,2)
  rho(i,j)=rhow+cf(i,j)*(sig-rhow)
  rhom(i,j)=rho(i,j)+cc(i,j)*(sig-rho(i,j))
end do
!$omp end parallel do
!
end subroutine cal_rho
!
!
!
subroutine r_iasc(fon,fname,ia,iend,jend,xllcorner,yllcorner,cellsize)
implicit none
integer :: j,iend,jend,fon,ncols,nrows,ia(1:iend,1:jend)
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
real(8) :: xllcorner,yllcorner,cellsize,a(1:iend,1:jend)
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
subroutine w_fasc_fmt(fon,fname,a,str,iend,jend,xllcorner,yllcorner,cellsize,nv)
implicit none
integer :: j,iend,jend,fon
real(8) :: a(1:iend,1:jend),xllcorner,yllcorner,cellsize,nv
character(len=100) :: fname,ciend,str
!
write(ciend,*)iend
open(fon,file=trim(adjustl(fname)))
write(fon,'(a,i5)')'ncols',iend
write(fon,'(a,i5)')'nrows',jend
write(fon,'(a,f15.3)')'xllcorner',xllcorner
write(fon,'(a,f15.3)')'yllcorner',yllcorner
write(fon,'(a,f10.3)')'cellsize ',cellsize
write(fon,'(a,f10.3)')'NODATA_value',nv
do j=jend,1,-1
  write(fon,trim(adjustl(str)))a(1:iend,j)
end do
close(fon)
!
end subroutine w_fasc_fmt
!
!
!
subroutine w_iasc_fmt(fon,fname,ia,str,iend,jend,xllcorner,yllcorner,cellsize,nv)
implicit none
integer :: j,iend,jend,fon,nv,ia(1:iend,1:jend)
real(8) :: xllcorner,yllcorner,cellsize
character(len=100) :: fname,ciend,str
!
write(ciend,*)iend
open(fon,file=trim(adjustl(fname)))
write(fon,'(a,i5)')'ncols',iend
write(fon,'(a,i5)')'nrows',jend
write(fon,'(a,f15.3)')'xllcorner',xllcorner
write(fon,'(a,f15.3)')'yllcorner',yllcorner
write(fon,'(a,f10.3)')'cellsize ',cellsize
write(fon,'(a,i5)')'NODATA_value',nv
do j=jend,1,-1
  write(fon,trim(adjustl(str)))ia(1:iend,j)
end do
close(fon)
!
end subroutine w_iasc_fmt
