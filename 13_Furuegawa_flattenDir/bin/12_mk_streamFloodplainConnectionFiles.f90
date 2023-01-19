implicit none
real(8), parameter :: PI=acos(-1.d0),deg2rad=PI/180.d0,rad2deg=180.d0/PI
!
integer :: n_ij_cv,n_ij_u,n_ij_v
integer, allocatable :: ij_cv(:,:),ij_u(:,:),ij_v(:,:),ibasin(:,:)
!
integer,allocatable :: ij_1d2d(:,:,:),i_to(:),i_from(:)
real(8),allocatable :: xy_1d2d(:,:)
real(8),allocatable :: f2d(:,:)
!
integer :: iend_1d,iofs,jofs,itr
integer :: i,j,iend,jend,itmp,k
real(8) :: x,y,xllcorner,yllcorner,cellsize,x1,x2,y1,y2
real(8) :: dist,dlx,dly,tmp,cs_c
character(len=100) :: fn,ctmp,ffmt
!
fn='targetArea_inDF.asc'
open(10,file=fn)
read(10,*)ctmp,iend
read(10,*)ctmp,jend
close(10)
write(*,*)iend,jend
allocate(ibasin(iend,jend))
call r_iasc(10,fn,ibasin,iend,jend,xllcorner,yllcorner,cellsize)
!
open(1,file='controlVolume_centerPoint_inDF.txt')
read(1,*)n_ij_cv
allocate(ij_cv(n_ij_cv,2))
read(1,*)
do k=1,n_ij_cv
  read(1,*)x,y,ij_cv(k,1:2)
end do
close(1)
!
open(1,file='fluxPoint_x_inDF.txt')
read(1,*)n_ij_u
allocate(ij_u(n_ij_u,2))
read(1,*)
do k=1,n_ij_u
  read(1,*)x,y,ij_u(k,1:2)
end do
close(1)
!
open(1,file='fluxPoint_y_inDF.txt')
read(1,*)n_ij_v
allocate(ij_v(n_ij_v,2))
read(1,*)
do k=1,n_ij_v
  read(1,*)x,y,ij_v(k,1:2)
end do
close(1)
!
open(10,file='streamConfiguration_inRR.txt')
read(10,*)iend_1d
allocate(ij_1d2d(iend_1d,2,1),xy_1d2d(iend_1d,2),i_to(iend_1d),i_from(iend_1d))
read(10,*)
do i=1,iend_1d
  read(10,*)itmp,i_to(i),i_from(i),itmp,itmp,tmp,tmp,tmp,dist,tmp,tmp,itmp,&
            xy_1d2d(i,1),xy_1d2d(i,2),ij_1d2d(i,1,1),ij_1d2d(i,2,1)
  if(i == 1)then
    cs_c=dist
  else if(dist < cs_c)then
    cs_c=dist
  end if
  !write(*,*)xy_1d2d(i,1),xy_1d2d(i,2),ij_1d2d(i,1),ij_1d2d(i,2),i_to(i)
end do
deallocate(ij_1d2d)
close(10)
!
itr=int(cs_c/cellsize)
allocate(ij_1d2d(iend_1d,2,itr))
ij_1d2d(:,:,:)=0
write(*,*)'# of div. =',itr
open(20,file='streamFloodplainConnection.txt')
write(20,'(i5,a)')itr,',# of div'
write(20,'(a)')'x,y,i2d,j2d,i1d,itr,dist'
do i=1,iend_1d
  !dist=(cellsize*dble(iend)*cellsize*dble(jend))**0.5d0
  do k=1,n_ij_cv
    x=xllcorner+cellsize*(dble(ij_cv(k,1)-1)+0.5d0)
    y=yllcorner+cellsize*(dble(ij_cv(k,2)-1)+0.5d0)
    dlx=xy_1d2d(i,1)-x
    dly=xy_1d2d(i,2)-y
    if(k==1)dist=(dlx*dlx+dly*dly)**0.5d0
    if((dlx*dlx+dly*dly)**0.5d0 <= dist)then
      ij_1d2d(i,1,1)=ij_cv(k,1)
      ij_1d2d(i,2,1)=ij_cv(k,2)
      dist=(dlx*dlx+dly*dly)**0.5d0
    end if
  end do
  if(ij_1d2d(i,1,1)==0 .or. ij_1d2d(i,2,1)==0)then
    write(*,*)'vol 1d -> 2d point missing @1D(i)',i,'<<<STOP>>>'
    stop
  end if
  !
  x1=xy_1d2d(i,1)!x1
  y1=xy_1d2d(i,2)!y1
  x2=xy_1d2d(i_to(i),1)!x2
  y2=xy_1d2d(i_to(i),2)!y2
  !
  if(i_to(i)==0)then
    x1=xy_1d2d(i_from(i),1)!x1
    y1=xy_1d2d(i_from(i),2)!y1
    x2=xy_1d2d(i,1)!x2
    y2=xy_1d2d(i,2)!y2
  end if
  !
  if(x2==x1)then
    iofs=0
  else
    iofs=nint(sign(1.d0,x2-x1))
  end if
  if(y2==y1)then
    jofs=0
  else
    jofs=nint(sign(1.d0,y2-y1))
  end if
  !
  do j=1,itr
    if(j==1)then
    else
      ij_1d2d(i,1,j)=ij_1d2d(i,1,j-1)+iofs
      ij_1d2d(i,2,j)=ij_1d2d(i,2,j-1)+jofs
    end if
    if(dist<cellsize/2.d0**0.5d0)then
      if(ibasin(ij_1d2d(i,1,j),ij_1d2d(i,2,j))==0)then
        write(*,*)'outof basin? @2D(i,j)',ij_1d2d(i,1,j),ij_1d2d(i,2,j)!,'<<<STOP>>>'
        !stop
      else
        write(20,'(2(f15.3,a),2(i5,a),(i10,a),(i5,a),f10.3)')&
          xllcorner+cellsize*(dble(ij_1d2d(i,1,j)-1)+0.5d0),',', &
          yllcorner+cellsize*(dble(ij_1d2d(i,2,j)-1)+0.5d0),',', &
          ij_1d2d(i,1,j),',',ij_1d2d(i,2,j),',',i,',',j,',',dist
      end if
    else
      !write(*,*)'dist=',dist,'<<<STOP>>>'
      !stop
    end if
  end do
end do
close(20)
!
allocate(f2d(iend,jend))
f2d=0.d0
fn='d.asc'
write(ctmp,*)maxval(f2d)
write(ffmt,*)len(trim(adjustl(ctmp)))+1
ffmt='(*(f'//trim(adjustl(ffmt))//'))'
ffmt='(*(f10.3))'
write(*,*)trim(adjustl(fn)),' ',ffmt
call w_fasc_fmt(20,fn,f2d,ffmt,iend,jend,xllcorner,yllcorner,cellsize,0.d0)
!
f2d=0.d0
fn='hs_ini.asc'
write(ctmp,*)maxval(f2d)
write(ffmt,*)len(trim(adjustl(ctmp)))+1
ffmt='(*(f'//trim(adjustl(ffmt))//'))'
ffmt='(*(f10.3))'
write(*,*)trim(adjustl(fn)),' ',ffmt
call w_fasc_fmt(20,fn,f2d,ffmt,iend,jend,xllcorner,yllcorner,cellsize,0.d0)
!
write(*,'(a)')'--- nomal end ---'
stop
end
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
