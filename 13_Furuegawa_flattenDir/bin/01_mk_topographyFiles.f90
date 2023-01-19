implicit none
integer,allocatable :: ibasin(:,:),iacc(:,:),idir(:,:),ibasin2(:,:),iacc2(:,:), &
                       i_2d(:,:),j_2d(:,:),id_2d(:,:)
integer,allocatable :: i_1d(:),j_1d(:),i_out(:),j_out(:)
integer :: iend,jend,i,j,k,i1,j1,i2,j2,itmp,itr,i0,j0,iout,jout
integer :: n_1d,n_out,n_ij_out,l
real(8),allocatable :: z_2d(:,:),facc(:,:)
real(8) :: xllcorner,yllcorner,cellsize,x_out,y_out
character(len=100) :: fname,ctmp,i_fmt,f_fmt
!
open(1,file='topographyConfiguration.txt')
read(1,*)iend
read(1,*)jend
allocate(ibasin(iend,jend),iacc(iend,jend),idir(iend,jend),ibasin2(iend,jend),iacc2(iend,jend), &
          i_2d(iend,jend),j_2d(iend,jend),id_2d(iend,jend), z_2d(iend,jend),facc(iend,jend))
read(1,*)
read(1,*)fname !dem
call r_fasc(10,fname,z_2d,iend,jend,xllcorner,yllcorner,cellsize)
fname='elevation_inRR.asc'
!write(ctmp,*)n_out
!write(i_fmt,*)len(trim(adjustl(ctmp)))+1
!i_fmt='(*(i'//trim(adjustl(i_fmt))//'))'
f_fmt='(*(f10.3)))'
!call w_iasc(20,fname,iriver(:,:),iend,jend,xllcorner,yllcorner,cellsize,0)
call w_fasc_fmt(20,fname,z_2d,f_fmt,iend,jend,xllcorner,yllcorner,cellsize,0.d0)
read(1,*)fname !dir
call r_iasc(10,fname,idir,iend,jend,xllcorner,yllcorner,cellsize)
idir=abs(idir)
!where(idir==0)idir=8
!where(idir==128)idir=0
!
ibasin=1
!
n_ij_out=1
allocate(i_out(n_ij_out),j_out(n_ij_out))
read(1,*)x_out
read(1,*)y_out
i_out=nint((x_out-xllcorner-0.5*cellsize)/cellsize)+1
j_out=nint((y_out-yllcorner-0.5*cellsize)/cellsize)+1
!
!! ---> find outlet
!n_ij_out=0
!do i=1,iend
!  do j=1,jend
!    if(1<=ibasin(i,j))then
!      call f_next(iend,jend,idir,i,j,i1,j1)
!      if(ibasin(i1,j1)==0)then
!        write(*,'(a,2i5)')'outlet(i,j)=',i,j
!        iout=i
!        jout=j
!        n_ij_out=n_ij_out+1
!      end if
!    end if
!  end do
!end do
!allocate(i_out(n_ij_out),j_out(n_ij_out))
!!
!n_ij_out=0
!do i=1,iend
!  do j=1,jend
!    if(1<=ibasin(i,j))then
!      call f_next(iend,jend,idir,i,j,i1,j1)
!      if(ibasin(i1,j1)==0)then
!        n_ij_out=n_ij_out+1
!        i_out(n_ij_out)=i
!        j_out(n_ij_out)=j
!      end if
!    end if
!  end do
!end do
!! <--- find outlet
!
i_2d(:,:)=0
j_2d(:,:)=0
id_2d(:,:)=0
ibasin2(:,:)=0
iacc2(:,:)=0
itr=1
do i=1,iend
  do j=1,jend
    i_2d(i,j)=i
    j_2d(i,j)=j
    id_2d(i,j)=(i-1)*jend+j
    if(itr/=(i-1)*jend+j)then
      write(*,*) itr,(i-1)*jend+j
      read(*,*)
    end if
    itr=itr+1
  end do
end do
!
n_1d=count(1<=ibasin(:,:))
allocate(i_1d(n_1d),j_1d(n_1d))
i_1d(1:n_1d)=pack(i_2d(:,:),(1<=ibasin(:,:)))
j_1d(1:n_1d)=pack(j_2d(:,:),(1<=ibasin(:,:)))
!write(*,'(a,i18)')' no. of data=',n_1d
!
do k=1,n_1d
  itmp=0
  i0=i_1d(k)
  j0=j_1d(k)
  i=i_1d(k)
  j=j_1d(k)
  itr=1
  i2=0
  j2=0
  do while(itmp==0)
    iacc2(i,j)=iacc2(i,j)+1
    call f_next(iend,jend,idir,i,j,i1,j1)
    do l=1,n_ij_out
      if(i==i_out(l) .and. j==j_out(l))then
        ibasin2(i0,j0)=1
        itmp=1
        !cycle
      end if
    end do
    if(i1<=1 .or. iend<=i1 .or. j1<=1 .or. jend<=j1)then
      itmp=1
    end if
    i2=i
    j2=j
    i=i1
    j=j1
    itr=itr+1
  end do
end do ! k
!
deallocate(i_1d,j_1d)
!
fname='targetArea_inRR.asc'
write(ctmp,*)n_out
write(i_fmt,*)len(trim(adjustl(ctmp)))+1
i_fmt='(*(i'//trim(adjustl(i_fmt))//'))'
!call w_iasc(20,fname,iriver(:,:),iend,jend,xllcorner,yllcorner,cellsize,0)
call w_iasc_fmt(20,fname,ibasin2,i_fmt,iend,jend,xllcorner,yllcorner,cellsize,0)
!
fname='flowAcc_inRR.asc'
write(ctmp,*)maxval(iacc2)
write(i_fmt,*)len(trim(adjustl(ctmp)))+1
i_fmt='(*(i'//trim(adjustl(i_fmt))//'))'
!call w_iasc(20,fname,iriver(:,:),iend,jend,xllcorner,yllcorner,cellsize,0)
call w_iasc_fmt(20,fname,iacc2,i_fmt,iend,jend,xllcorner,yllcorner,cellsize,0)
!
fname='flowDir_inRR.asc'
write(ctmp,*)maxval(idir)
write(i_fmt,*)len(trim(adjustl(ctmp)))+1
i_fmt='(*(i'//trim(adjustl(i_fmt))//'))'
!call w_iasc(20,fname,iriver(:,:),iend,jend,xllcorner,yllcorner,cellsize,0)
call w_iasc_fmt(20,fname,idir,i_fmt,iend,jend,xllcorner,yllcorner,cellsize,0)
!
fname='volcanicAsh_inRR.asc'
write(ctmp,*)maxval(idir)
write(f_fmt,*)len(trim(adjustl(ctmp)))+1
f_fmt='(*(f10.3)))'
!call w_iasc(20,fname,iriver(:,:),iend,jend,xllcorner,yllcorner,cellsize,0)
call w_fasc_fmt(20,fname,dble(idir)*0.d0,f_fmt,iend,jend,xllcorner,yllcorner,cellsize,0.d0)
!
write(*,'(a)')'--- nomal end ---'
stop
end
!
!
!
subroutine b_up(iend,jend,iriver,idir,i,j,itmp)
implicit none
integer :: i,j,iend,jend,itmp
integer :: iriver(iend,jend),idir(iend,jend)
!
! nw n ne ! 3  2  1
!  \ | /  !  \ | /
! w -+- e ! 4 -+- 8
!  / | \　　!  / | \
! sw s se　! 5  6 　7　　　　
!if(iriver(i,  j+1)==1)itmp=itmp+1 ! n
!if(iriver(i+1,j+1)==1)itmp=itmp+1 ! ne
!if(iriver(i+1,  j)==1)itmp=itmp+1 ! e
!if(iriver(i+1,j-1)==1)itmp=itmp+1 ! se
!if(iriver(i,  j-1)==1)itmp=itmp+1 ! s
!if(iriver(i-1,j-1)==1)itmp=itmp+1 ! sw
!if(iriver(i-1,  j)==1)itmp=itmp+1 ! w
!if(iriver(i-1,j+1)==1)itmp=itmp+1 ! nw
if(idir(i-1,j-1)==1 .and. iriver(i-1,j-1)==1)itmp=itmp+1 ! 5 -> 1
if(idir(i  ,j-1)==2 .and. iriver(i  ,j-1)==1)itmp=itmp+1 ! 6 -> 2
if(idir(i+1,j-1)==3 .and. iriver(i+1,j-1)==1)itmp=itmp+1 ! 7 -> 3
if(idir(i+1,j  )==4 .and. iriver(i+1,j  )==1)itmp=itmp+1 ! 8 -> 4
if(idir(i+1,j+1)==5 .and. iriver(i+1,j+1)==1)itmp=itmp+1 ! 1 -> 5
if(idir(i  ,j+1)==6 .and. iriver(i  ,j+1)==1)itmp=itmp+1 ! 2 -> 6
if(idir(i-1,j+1)==7 .and. iriver(i-1,j+1)==1)itmp=itmp+1 ! 3 -> 7
if(idir(i-1,j  )==8 .and. iriver(i-1,j  )==1)itmp=itmp+1 ! 4 -> 8
!
end subroutine b_up
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
!else if(idir(i,j)<=-1 .or. iriver(i,j)==0)then
!  write(*,'(a,4i5)')'river end @ i,j,dir,iriv=',i,j,idir(i,j),iriver(i,j)
!  itmp=1
else if(idir(i,j)==0)then
  i1=-1
  j1=-1
else
  write(*,'(a,4i5)')'river course? @ i,j,dir,iriv=',i,j,idir(i,j) !,iriver(i,j)
  write(*,*)'Don''t come here?'
  write(*,*)'iend,jend',iend,jend
  stop
end if
!
end subroutine f_next
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
subroutine w_iasc(fon,fname,ia,iend,jend,xllcorner,yllcorner,cellsize,nv)
implicit none
integer :: j,iend,jend,fon,nv,ia(1:iend,1:jend)
real(8) :: xllcorner,yllcorner,cellsize
character(len=100) :: fname,ciend
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
  write(fon,'(*(i2))')ia(1:iend,j)
end do
close(fon)
!
end subroutine w_iasc
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
subroutine w_fasc(fon,fname,a,iend,jend,xllcorner,yllcorner,cellsize,nv)
implicit none
integer :: j,iend,jend,fon,nv
real(8) :: xllcorner,yllcorner,cellsize,a(1:iend,1:jend)
character(len=100) :: fname,ciend
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
  write(fon,'(*(f10.3))')a(1:iend,j)
end do
close(fon)
!
end subroutine w_fasc
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
