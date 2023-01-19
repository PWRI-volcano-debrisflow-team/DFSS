implicit none
integer, allocatable :: ibasin(:,:)
integer, allocatable :: i_2d(:,:),j_2d(:,:)
integer, allocatable :: iqx_bd_basin(:,:),iqy_bd_basin(:,:)
! scalar 2d -> 1d
logical, allocatable :: m_cv(:,:)
integer :: n_ij_cv
integer, allocatable :: ij_cv(:,:)
! vector 2d -> 1d
logical, allocatable :: m_u(:,:),m_v(:,:)
integer :: n_ij_u,n_ij_v
integer, allocatable :: ij_u(:,:),ij_v(:,:)
! basin boundary
logical, allocatable :: m_u_w(:,:),m_u_e(:,:),m_v_n(:,:),m_v_s(:,:)
integer :: n_ij_u_w,n_ij_u_e,n_ij_v_n,n_ij_v_s, &
           n_ij_u_we,n_ij_v_sn
integer, allocatable :: ij_u_w(:,:),ij_u_e(:,:), ij_v_n(:,:),ij_v_s(:,:), &
                        ij_u_we(:,:),ij_v_sn(:,:)
integer :: i,iend,j,jend
real(8) :: xllcorner,yllcorner,cellsize,dx,dy
character(len=100) :: fname,head,fn_bsn,fn_dir,fn_dem
!
real(8), parameter :: PI=acos(-1.d0),D2R=PI/180.d0,R2D=180.d0/PI
integer :: k,i1,j1
real(8) :: dl,deg,tele,dtr,dep
integer, allocatable :: idir(:,:)
real(8), allocatable :: z(:,:),d(:,:),grad(:,:),hsc(:,:),str(:,:)
character(len=100) :: ctmp,i_fmt,f_fmt
!
open(1,file='topographyConfiguration_inDF.txt')
read(1,*)iend
read(1,*)jend
read(1,*)
read(1,*)fn_bsn
read(1,*)fn_dem
read(1,*)fn_dir
close(1)
!
allocate(ibasin(iend,jend),z(iend,jend),idir(iend,jend))
allocate(i_2d(iend+1,jend+1),j_2d(iend+1,jend+1))
allocate(iqx_bd_basin(iend+1,jend+1),iqy_bd_basin(iend+1,jend+1))
!
call r_iasc(10,fn_bsn,ibasin,iend,jend,xllcorner,yllcorner,cellsize)
fname='targetArea_inDF.asc'
write(ctmp,*)iend
write(f_fmt,*)len(trim(adjustl(ctmp)))+1
f_fmt='(*(i'//trim(adjustl(f_fmt))//'))'
!call w_iasc(20,fname,iriver(:,:),iend,jend,xllcorner,yllcorner,cellsize,0)
call w_iasc_fmt(20,fname,ibasin,f_fmt,iend,jend,xllcorner,yllcorner,cellsize,0)
!
call r_fasc(10,fn_dem,z,iend,jend,xllcorner,yllcorner,cellsize)
fname='elevation_inDF.asc'
f_fmt='(*(f10.3)))'
call w_fasc_fmt(20,fname,z,f_fmt,iend,jend,xllcorner,yllcorner,cellsize,-999.d0)
!
call r_iasc(10,fn_dir,idir,iend,jend,xllcorner,yllcorner,cellsize)
fname='flowDir_inDF.asc'
write(ctmp,*)maxval(idir)
write(i_fmt,*)len(trim(adjustl(ctmp)))+1+1
i_fmt='(*(i'//trim(adjustl(i_fmt))//'))'
call w_iasc_fmt(20,fname,idir,i_fmt,iend,jend,xllcorner,yllcorner,cellsize,0)
!
dx=cellsize
dy=cellsize
!
do i=1,iend+1
  do j=1,jend+1
    i_2d(i,j)=i
    j_2d(i,j)=j
  end do
end do
! ---> set control volume
allocate(m_cv(iend,jend))
m_cv(1:iend,1:jend)=(ibasin(1:iend,1:jend)/=0)
n_ij_cv=count(m_cv)
allocate(ij_cv(n_ij_cv,2))
ij_cv(:,1)=pack(i_2d(1:iend,1:jend),m_cv)
ij_cv(:,2)=pack(j_2d(1:iend,1:jend),m_cv)
! <--- set control volume
!
! ---> set boundary flag
iqx_bd_basin(1:iend+1,1:jend+1)=-1
iqy_bd_basin(1:iend+1,1:jend+1)=-1
call set_flux_boundary_basin(iend,jend,ibasin(:,:), &
      iqx_bd_basin(1:iend+1,1:jend+1),iqy_bd_basin(1:iend+1,1:jend+1))
! boundary flags
! 1: flux@slope qx,qy
! 50: north basin qy
! 51: east  basin qx
! 52: south basin qy
! 53: west  basin qx
! <--- set boundary flag

! ---> set flux
allocate(m_u(iend+1,jend+1), m_v(iend+1,jend+1))
!m_u(1:iend+1,1:jend+1)=(iqx_bd_basin(1:iend+1,1:jend+1)==1) ! mask
!m_v(1:iend+1,1:jend+1)=(iqy_bd_basin(1:iend+1,1:jend+1)==1) ! mask
m_u(1:iend+1,1:jend+1)=(iqx_bd_basin(1:iend+1,1:jend+1)==1  .or. &
                        iqx_bd_basin(1:iend+1,1:jend+1)==51 .or. &
                        iqx_bd_basin(1:iend+1,1:jend+1)==53) ! mask
m_v(1:iend+1,1:jend+1)=(iqy_bd_basin(1:iend+1,1:jend+1)==1  .or. &
                        iqy_bd_basin(1:iend+1,1:jend+1)==50 .or. &
                        iqy_bd_basin(1:iend+1,1:jend+1)==52) ! mask
n_ij_u=count(m_u)
n_ij_v=count(m_v)
allocate(ij_u(n_ij_u,2),ij_v(n_ij_v,2))
ij_u(:,1)=pack(i_2d,m_u)
ij_u(:,2)=pack(j_2d,m_u)
ij_v(:,1)=pack(i_2d,m_v)
ij_v(:,2)=pack(j_2d,m_v)
! <---  set flux
!
! ---> set boundary
allocate(m_u_w(iend+1,jend+1),m_u_e(iend+1,jend+1), &
         m_v_n(iend+1,jend+1),m_v_s(iend+1,jend+1))
m_u_w(1:iend+1,1:jend+1)=(iqx_bd_basin(1:iend+1,1:jend+1)==53)! mask
m_u_e(1:iend+1,1:jend+1)=(iqx_bd_basin(1:iend+1,1:jend+1)==51)! mask
m_v_s(1:iend+1,1:jend+1)=(iqy_bd_basin(1:iend+1,1:jend+1)==52)! mask
m_v_n(1:iend+1,1:jend+1)=(iqy_bd_basin(1:iend+1,1:jend+1)==50)! mask
n_ij_u_w=count(m_u_w)
n_ij_u_e=count(m_u_e)
n_ij_v_s=count(m_v_s)
n_ij_v_n=count(m_v_n)
allocate(ij_u_w(n_ij_u_w,2), &
         ij_u_e(n_ij_u_e,2), &
         ij_v_s(n_ij_v_s,2), &
         ij_v_n(n_ij_v_n,2))
ij_u_w(:,1)=pack(i_2d,m_u_w)
ij_u_w(:,2)=pack(j_2d,m_u_w)
ij_u_e(:,1)=pack(i_2d,m_u_e)
ij_u_e(:,2)=pack(j_2d,m_u_e)
ij_v_s(:,1)=pack(i_2d,m_v_s)
ij_v_s(:,2)=pack(j_2d,m_v_s)
ij_v_n(:,1)=pack(i_2d,m_v_n)
ij_v_n(:,2)=pack(j_2d,m_v_n)
n_ij_u_we=n_ij_u_w+n_ij_u_e
n_ij_v_sn=n_ij_v_s+n_ij_v_n
!
allocate(ij_u_we(n_ij_u_we,8),ij_v_sn(n_ij_v_sn,8))
!
do i=1,n_ij_u_w
  !west  h(i  ,  j),(wl(i+1,  j)-wl(i  ,  j))/dx (0 , 0),( 1, 0),( 0, 0),flux<0
  ij_u_we(i,1:2)=ij_u_w(i,1:2)
  ij_u_we(i,3:8)=(/0,0,1,0,0,0/)
end do
do i=1,n_ij_u_e
  !east  h(i-1,  j),(wl(i-1,  j)-wl(i-2,  j))/dx (-1, 0),(-1, 0),(-2, 0),flux>0
  ij_u_we(n_ij_u_w+i,1:2)=ij_u_e(i,1:2)
  ij_u_we(n_ij_u_w+i,3:8)=(/-1,0,-1,0,-2,0/)
end do
do i=1,n_ij_v_s
  !south h(i  ,  j),(wl(i  ,j+1)-wl(i  ,  j))/dy (0 , 0),( 0, 1),( 0, 0),flux<0
  ij_v_sn(i,1:2)=ij_v_s(i,1:2)
  ij_v_sn(i,3:8)=(/0,0,0,1,0,0/)
end do
do i=1,n_ij_v_n
  !north h(i  ,j-1),(wl(i  ,j-1)-wl(i  ,j-2))/dy (0 ,-1),( 0,-1),( 0,-2),flux>0
  ij_v_sn(n_ij_v_s+i,1:2)=ij_v_n(i,1:2)
  ij_v_sn(n_ij_v_s+i,3:8)=(/0,-1,0,-1,0,-2/)
end do
! <--- set boundary
!
! ---> output
fname=trim(adjustl('controlVolume_centerPoint_inDF.txt'))
call point2(20,fname,xllcorner+dx*0.5,yllcorner+dy*0.5,dx,n_ij_cv,ij_cv(:,1),ij_cv(:,2))
fname=trim(adjustl('fluxPoint_x_inDF.txt'))
call point2(20,fname,xllcorner,yllcorner+dy*0.5,dx,n_ij_u,ij_u(:,1),ij_u(:,2))
fname=trim(adjustl('fluxPoint_y_inDF.txt'))
call point2(20,fname,xllcorner+dx*0.5,yllcorner,dx,n_ij_v,ij_v(:,1),ij_v(:,2))
!
fname=trim(adjustl('fluxPoint_xBoudary_inDF.txt'))
head='x,y,cvx,cvy,f1x,f1y,f2x,f2y,f3x,f3y,cs'
call point3(20,fname,xllcorner,yllcorner+dy*0.5,cellsize, &
            n_ij_u_w+n_ij_u_e,ij_u_we(:,1),ij_u_we(:,2),ij_u_we(:,1:),8,head)
fname=trim(adjustl('fluxPoint_yBoudary_inDF.txt'))
head='x,y,cvx,cvy,f1x,f1y,f2x,f2y,f3x,f3y,cs'
call point3(20,fname,xllcorner+dx*0.5,yllcorner,cellsize, &
            n_ij_v_s+n_ij_v_n,ij_v_sn(:,1),ij_v_sn(:,2),ij_v_sn(:,1:),8,head)
! <--- output
!
open(10,file='fluxPoint_wBoudary_inDF.txt')
write(10,*)n_ij_u_w
write(10,'(a)')'x,y,i,j'
do k=1,n_ij_u_w
  i=ij_u_w(k,1)
  j=ij_u_w(k,2)
  write(10,'(2(f15.3,:,","),2(i5,:,","))') &
  xllcorner+cellsize*float(i-1),yllcorner+dy*0.5+cellsize*float(j-1),i,j
end do
close(10)
!
open(10,file='fluxPoint_eBoudary_inDF.txt')
write(10,*)n_ij_u_e
write(10,'(a)')'x,y,i,j'
do k=1,n_ij_u_e
  i=ij_u_e(k,1)
  j=ij_u_e(k,2)
  write(10,'(2(f15.3,:,","),2(i5,:,","))') &
  xllcorner+cellsize*float(i-1),yllcorner+dy*0.5+cellsize*float(j-1),i,j
end do
close(10)
!
open(10,file='fluxPoint_sBoudary_inDF.txt')
write(10,*)n_ij_v_s
write(10,'(a)')'x,y,i,j'
do k=1,n_ij_v_s
  i=ij_v_s(k,1)
  j=ij_v_s(k,2)
  write(10,'(2(f15.3,:,","),2(i5,:,","))') &
  xllcorner+dx*0.5+cellsize*float(i-1),yllcorner+cellsize*float(j-1),i,j
end do
close(10)
!
open(10,file='fluxPoint_nBoudary_inDF.txt')
write(10,*)n_ij_v_n
write(10,'(a)')'x,y,i,j'
do k=1,n_ij_v_n
  i=ij_v_n(k,1)
  j=ij_v_n(k,2)
  write(10,'(2(f15.3,:,","),2(i5,:,","))') &
  xllcorner+dx*0.5+cellsize*float(i-1),yllcorner+cellsize*float(j-1),i,j
end do
close(10)
!
do i=1,n_ij_u_e
  !east  h(i-1,  j),(wl(i-1,  j)-wl(i-2,  j))/dx (-1, 0),(-1, 0),(-2, 0),flux>0
  ij_u_we(n_ij_u_w+i,1:2)=ij_u_e(i,1:2)
  ij_u_we(n_ij_u_w+i,3:8)=(/-1,0,-1,0,-2,0/)
end do
do i=1,n_ij_v_s
  !south h(i  ,  j),(wl(i  ,j+1)-wl(i  ,  j))/dy (0 , 0),( 0, 1),( 0, 0),flux<0
  ij_v_sn(i,1:2)=ij_v_s(i,1:2)
  ij_v_sn(i,3:8)=(/0,0,0,1,0,0/)
end do
do i=1,n_ij_v_n
  !north h(i  ,j-1),(wl(i  ,j-1)-wl(i  ,j-2))/dy (0 ,-1),( 0,-1),( 0,-2),flux>0
  ij_v_sn(n_ij_v_s+i,1:2)=ij_v_n(i,1:2)
  ij_v_sn(n_ij_v_s+i,3:8)=(/0,-1,0,-1,0,-2/)
end do
!
write(*,'(a)')'--- nomal end ---'
stop
end
!
!
!
subroutine set_flux_boundary_basin(iend,jend,ibasin,iqx_bd_rs,iqy_bd_rs)
implicit none
integer :: i,j,iend,jend
integer :: ibasin(1:iend,1:jend),iqx_bd_rs(1:iend+1,1:jend+1),iqy_bd_rs(1:iend+1,1:jend+1)
!
! boundary flags
! 50: north basin qy
! 51: east  basin qx
! 52: south basin qy
! 53: west  basin qx
do i=2,iend
  do j=2,jend
    if(ibasin(i-1,j)==0 .and. ibasin(i,j)/=0)then
      iqx_bd_rs(i,j)=53 !west
    else if(ibasin(i-1,j)/=0 .and. ibasin(i,j)==0)then
      iqx_bd_rs(i,j)=51 !east
    else if(ibasin(i-1,j)/=0 .and. ibasin(i,j)/=0)then
      iqx_bd_rs(i,j)=1
    end if
    if(ibasin(i,j-1)==0 .and. ibasin(i,j)/=0)then
      iqy_bd_rs(i,j)=52 !south
    else if(ibasin(i,j-1)/=0 .and. ibasin(i,j)==0)then
      iqy_bd_rs(i,j)=50 !north
    else if(ibasin(i,j-1)/=0 .and. ibasin(i,j)/=0)then
      iqy_bd_rs(i,j)=1
    end if
  end do
end do
!
do i=2,iend
  j=1
  if(ibasin(i-1,j)==0 .and. ibasin(i,j)/=0)then
    iqx_bd_rs(i,j)=53 !west
  else if(ibasin(i-1,j)/=0 .and. ibasin(i,j)==0)then
    iqx_bd_rs(i,j)=51 !east
  end if
end do
!
do j=2,jend
i=1
  if(ibasin(i,j-1)==0 .and. ibasin(i,j)/=0)then
    iqy_bd_rs(i,j)=52 !south
  else if(ibasin(i,j-1)/=0 .and. ibasin(i,j)==0)then
    iqy_bd_rs(i,j)=50 !north
  end if
end do
! ---> qx boundary
i=1
do j=1,jend
  if(ibasin(i,j)/=0)then
    iqx_bd_rs(i,j)=53 !west
  end if
end do
!
i=iend+1
do j=1,jend
  if(ibasin(i-1,j)/=0)then
    iqx_bd_rs(i,j)=51 !east
  end if
end do
! <--- qx boundary
!
! ---> qy boundary
j=1
do i=1,iend
  if(ibasin(i,j)/=0)then
    iqy_bd_rs(i,j)=52 !south
  end if
end do
!
j=jend+1
do i=1,iend
  if(ibasin(i,j-1)/=0)then
    iqy_bd_rs(i,j)=50 !north
  end if
end do
! <--- qy boundary
!
end subroutine set_flux_boundary_basin
!
!
!
subroutine w_fasc(fon,fname,cfmt,a,iend,jend,xllcorner,yllcorner,cellsize,nv)
implicit none
integer :: j,iend,jend,fon
real(8) :: xllcorner,yllcorner,cellsize,nv,a(1:iend,1:jend)
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
integer :: j,iend,jend,fon,nv,ia(1:iend,1:jend)
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
subroutine point(fon,fname,a,iend,jend,xllcorner,yllcorner,cellsize)
implicit none
integer :: i,j,iend,jend,fon,a(1:iend,1:jend),iflg
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
write(fon,*)n_1d
write(fon,'(a)')'x,y,i,j,cs'
do ni=1,n_1d
  i=i_1d(ni)
  j=j_1d(ni)
  write(fon,'(2(f15.3,a),2(i5,a),f10.3)')&
  xllcorner+cellsize*float(i-1),',',yllcorner+cellsize*float(j-1),',',i,',',j,',',cellsize
end do
close(fon)
!
end subroutine point2
!
!
!
subroutine point3(fon,fname,xllcorner,yllcorner,cellsize,n_1d,i_1d,j_1d,flg,n_flg,head)
implicit none
integer :: ni,n_1d,n_flg,i_1d(n_1d),j_1d(n_1d),flg(n_1d,n_flg),i,j,fon
real(8) :: xllcorner,yllcorner,cellsize
character(len=100) :: head,fname,cn_flg
!
open(fon,file=trim(adjustl(fname)))
write(fon,*)n_1d
write(fon,'(a)')head
write(cn_flg,*)n_flg
do ni=1,n_1d
  i=i_1d(ni)
  j=j_1d(ni)
  write(fon,'(2(f15.3,:,","),'//trim(adjustl(cn_flg))//'(i5,:,","),f10.3)') &
  xllcorner+cellsize*float(i-1),yllcorner+cellsize*float(j-1),flg(ni,:),cellsize
end do
close(fon)
!
end subroutine point3
!
!
!
subroutine f_next(iend,jend,idir,i,j,i1,j1)
implicit none
integer :: i,j,iend,jend,i1,j1,idir(iend,jend)
!
! nw n ne ! 3  2  1
!  \ | /  !  \ | /
! w -+- e ! 4 -+- 8
!  / | \　　!  / | \
! sw s se　! 5  6 　7　　　　
if(     idir(i,j)==1)then
  i1=i+1
  j1=j+1
else if(idir(i,j)==2)then
  i1=i
  j1=j+1
else if(idir(i,j)==3)then
  i1=i-1
  j1=j+1
else if(idir(i,j)==4)then
  i1=i-1
  j1=j
else if(idir(i,j)==5)then
  i1=i-1
  j1=j-1
else if(idir(i,j)==6)then
  i1=i
  j1=j-1
else if(idir(i,j)==7)then
  i1=i+1
  j1=j-1
else if(idir(i,j)==8)then
  i1=i+1
  j1=j
else if(i1<1 .or. i1>iend .or. j1<1 .or. j1>jend)then
  i1=-1
  j1=-1
else
  write(*,'(a,4i5)')'river course? @ i,j,dir,iriv=',i,j,idir(i,j)
  write(*,*)'Don''t come here?'
  stop
end if
!
end subroutine f_next
!
!
!
subroutine cal_stability_hsc(n_1d,ij_1d,iend,jend,d,grad,hsc)
implicit none
real(8), parameter :: PI=acos(-1.d0),D2R=PI/180.d0,R2D=180.d0/PI
real(8), parameter :: g=9.8d0,sig=2650.d0,rho=1000.d0,lamda=0.4d0,pw=0.1d0, &
                      c=3000.d0, tanp=tan(35.d0*D2R)
integer :: ni,n_1d,ij_1d(n_1d,2),i,j,iend,jend
real(8) :: d(iend,jend),grad(iend,jend),c2,hsc(iend,jend),hsc0,csta,tant,cost
!
csta=1.d0-lamda
do ni=1,n_1d
  i=ij_1d(ni,1)
  j=ij_1d(ni,2)
  tant=grad(i,j)*D2R
  cost=cos(atan(tant))
  c2=c/(rho*g*d(i,j)*cost*tanp)
  hsc0=((1.d0-tant/tanp)*(csta*sig/rho+pw)+c2) / &
       ((1.d0-tant/tanp)*(csta+pw)+tant/tanp)
  hsc(i,j)=hsc0/d(i,j)
end do
!
end subroutine cal_stability_hsc
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
  write(fon,str)ia(1:iend,j)
end do
close(fon)
!
end subroutine w_iasc_fmt
!
!
!
subroutine w_fasc_fmt(fon,fname,a,str,iend,jend,xllcorner,yllcorner,cellsize,nv)
implicit none
integer :: j,iend,jend,fon
real(8) :: xllcorner,yllcorner,cellsize,nv,a(1:iend,1:jend)
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
  write(fon,str)a(1:iend,j)
end do
close(fon)
!
end subroutine w_fasc_fmt






