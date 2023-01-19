implicit none
real(8), parameter :: PI=acos(-1.d0),deg2rad=PI/180.d0,rad2deg=180.d0/PI
integer,parameter :: upcell=1 ! start point of river
integer,allocatable :: ibasin(:,:),idrain(:,:),iriver(:,:),istream_order(:,:),iacc(:,:), &
                       i_2d(:,:),j_2d(:,:),i1_2d(:,:),j1_2d(:,:),i2_2d(:,:),j2_2d(:,:)
integer :: iend,jend,i,j,k,i1,j1,i2,j2,itmp,iacc_thresh,nirs,itr,iso_max
integer,allocatable :: i_1d(:),j_1d(:),in_1d(:),jn_1d(:),ip_1d(:),jp_1d(:),idrain_1d(:),iso_1d(:),iacc_1d(:), &
                       iriver_sta_i(:),iriver_sta_j(:)
integer :: n_1d,n_inf,iflg_in0,iflg_in1
real(8),allocatable :: x_2d(:,:),y_2d(:,:),z_2d(:,:)!,acc(iend,jend)
real(8),allocatable :: x_1d(:),y_1d(:),z_1d(:)!,acc_1d(:)
real(8) :: x,y,xllcorner,yllcorner,cellsize,ba
character(len=100) :: fname,ctmp,ftmp
!
integer,allocatable :: id_from(:,:),id_to(:),idx(:)
integer :: iexit
logical :: mask
real(8),allocatable :: b_1d(:),d_1d(:),dx(:)
real(8) :: dx0,dy0,x_vec,y_vec,deg,dl,grad,wm,wp,wmin,dm,dp,dmin,mn,dsoil
!
open(1,file='topographyConfiguration.txt')
read(1,*)iend
read(1,*)jend
allocate(ibasin(iend,jend),idrain(iend,jend),iriver(iend,jend),istream_order(iend,jend),iacc(iend,jend), &
         i_2d(iend,jend),j_2d(iend,jend),i1_2d(iend,jend),j1_2d(iend,jend),i2_2d(iend,jend),j2_2d(iend,jend), &
         x_2d(iend,jend),y_2d(iend,jend),z_2d(iend,jend))
read(1,*)
read(1,*)
!
fname='elevation_inRR.asc'
call r_fasc(10,fname,z_2d,iend,jend,xllcorner,yllcorner,cellsize)
!
read(1,*)
read(1,*)
read(1,*)
read(1,*)
read(1,*)
read(1,*)iacc_thresh
!
fname='targetArea_inRR.asc'
call r_iasc(10,fname,ibasin,iend,jend,xllcorner,yllcorner,cellsize)
!
fname='flowAcc_inRR.asc'
call r_iasc(10,fname,iacc,iend,jend,xllcorner,yllcorner,cellsize)
where(iacc(1:iend,1:jend)>iacc_thresh .and. ibasin(1:iend,1:jend)/=0)
  iriver(1:iend,1:jend)=1
elsewhere
  iriver(1:iend,1:jend)=0
end where
!
fname='flowDir_inRR.asc'
call r_iasc(10,fname,idrain,iend,jend,xllcorner,yllcorner,cellsize)
!
! --> set start of river
! nw n ne ! 3  2  1
!  \ | /  !  \ | /
! w -+- e ! 4 -+- 8
!  / | \  !  / | \
! sw s se ! 5  6  7
nirs=0
open(20,file='river_bounds.txt')
write(20,*)'x,y,i,j'
do i=1,iend
  do j=1,jend
    if(iriver(i,j)==1)then
      !
      iflg_in1=0
      if(idrain(i,  j+1)==6 .and. iriver(i,  j+1)==1)iflg_in1=iflg_in1+1 ! n
      if(idrain(i+1,j+1)==5 .and. iriver(i+1,j+1)==1)iflg_in1=iflg_in1+1 ! ne
      if(idrain(i+1,  j)==4 .and. iriver(i+1,  j)==1)iflg_in1=iflg_in1+1 ! e
      if(idrain(i+1,j-1)==3 .and. iriver(i+1,j-1)==1)iflg_in1=iflg_in1+1 ! se
      if(idrain(i,  j-1)==2 .and. iriver(i,  j-1)==1)iflg_in1=iflg_in1+1 ! s
      if(idrain(i-1,j-1)==1 .and. iriver(i-1,j-1)==1)iflg_in1=iflg_in1+1 ! sw
      if(idrain(i-1,  j)==8 .and. iriver(i-1,  j)==1)iflg_in1=iflg_in1+1 ! w
      if(idrain(i-1,j+1)==7 .and. iriver(i-1,j+1)==1)iflg_in1=iflg_in1+1 ! nw
      !
      iflg_in0=0
      if(idrain(i,  j+1)==6 .and. iriver(i,  j+1)==0)iflg_in0=iflg_in0+1 ! n
      if(idrain(i+1,j+1)==5 .and. iriver(i+1,j+1)==0)iflg_in0=iflg_in0+1 ! ne
      if(idrain(i+1,  j)==4 .and. iriver(i+1,  j)==0)iflg_in0=iflg_in0+1 ! e
      if(idrain(i+1,j-1)==3 .and. iriver(i+1,j-1)==0)iflg_in0=iflg_in0+1 ! se
      if(idrain(i,  j-1)==2 .and. iriver(i,  j-1)==0)iflg_in0=iflg_in0+1 ! s
      if(idrain(i-1,j-1)==1 .and. iriver(i-1,j-1)==0)iflg_in0=iflg_in0+1 ! sw
      if(idrain(i-1,  j)==8 .and. iriver(i-1,  j)==0)iflg_in0=iflg_in0+1 ! w
      if(idrain(i-1,j+1)==7 .and. iriver(i-1,j+1)==0)iflg_in0=iflg_in0+1 ! nw
      !
      itmp=0
      if(iriver(i,  j+1)==1)itmp=itmp+1 ! n
      if(iriver(i+1,j+1)==1)itmp=itmp+1 ! ne
      if(iriver(i+1,  j)==1)itmp=itmp+1 ! e
      if(iriver(i+1,j-1)==1)itmp=itmp+1 ! se
      if(iriver(i,  j-1)==1)itmp=itmp+1 ! s
      if(iriver(i-1,j-1)==1)itmp=itmp+1 ! sw
      if(iriver(i-1,  j)==1)itmp=itmp+1 ! w
      if(iriver(i-1,j+1)==1)itmp=itmp+1 ! nw
      !
      if(iflg_in0>=upcell .and. iflg_in1==0)then
      !if(itmp==1  .or. iflg_in==0)then
        call f_next(iend,jend,idrain,i,j,i1,j1)
        if (iriver(i1,j1)==1)then
          nirs=nirs+1
         ! write(*,*)'inlet(i,j)=',i,j
          write(20,'(2(f15.3,a),3(i5,a))')xllcorner+cellsize*(0.5+i-1),',', &
          yllcorner+cellsize*(0.5+j-1),',',i,',',j
        else
         ! write(*,*)'outlet',i,j,'Ent.'
          !read(*,*)
        end if
      end if
    end if
  end do
end do
close(20)
!
allocate(iriver_sta_i(nirs),iriver_sta_j(nirs))
!write(*,*)'Input OK, Ent.','No. point =',nirs
! <-- set start of river
!
istream_order(:,:)=0
x_2d(:,:)=-1
y_2d(:,:)=-1
i_2d(:,:)=0
j_2d(:,:)=0
i1_2d(:,:)=0
j1_2d(:,:)=0
i2_2d(:,:)=0
j2_2d(:,:)=0
open(20,file='river_bounds.txt')
read(20,*)
do k=1,nirs
  read(20,*)ftmp,ftmp,iriver_sta_i(k),iriver_sta_j(k) !,nirs_dirflag(k)
  i=iriver_sta_i(k)
  j=iriver_sta_j(k)
  !write(ctmp,'(i3.3)')k
  !write(ctmp,'(a,i3.3,a,i3.3)')'i',i,'_j',j
  !fname='rivers/river_'//trim(adjustl(ctmp))//'.csv'
  !open(100,file=fname)
  !write(100,'(a)')'x,y,z,i,j,i-1,j-1,i+1,j+1'
  !write(*,*)iriver_sta_i(k),iriver_sta_j(k) !,nirs_dirflag(k)
  itmp=0
  !i=iriver_sta_i(k)
  !j=iriver_sta_j(k)
  istream_order(i,j)=1
  itr=1
  i2=0
  j2=0
  do while(itmp==0)
    call f_next(iend,jend,idrain,i,j,i1,j1)
    x=xllcorner+cellsize*(0.5+i-1)
    y=yllcorner+cellsize*(0.5+j-1)
    x_2d(i,j)=x
    y_2d(i,j)=y
    i_2d(i,j)=i
    j_2d(i,j)=j
    i1_2d(i,j)=i1
    j1_2d(i,j)=j1
    i2_2d(i,j)=i2
    j2_2d(i,j)=j2
    !if(iriver(i1,j1)==0)then
    !  write(100,'(3(f15.3,a),6(i5,a))')x,',',y,',',z_2d(i,j),',',i,',',j,',',i2,',',j2,',',0,',',0
    !else
    !  write(100,'(3(f15.3,a),6(i5,a))')x,',',y,',',z_2d(i,j),',',i,',',j,',',i2,',',j2,',',i1,',',j1
    !end if
    if(i1==-1 .or. j1==-1 .or. ibasin(i1,j1)==0)then
      exit
    end if
    !
    iso_max=0
    n_inf=0
    if(idrain(i1+1,j1+1)==5 .and. idrain(i,j)/=5) iso_max=max(iso_max,istream_order(i1+1,j1+1))
    if(idrain(i1  ,j1+1)==6 .and. idrain(i,j)/=6) iso_max=max(iso_max,istream_order(i1  ,j1+1))
    if(idrain(i1-1,j1+1)==7 .and. idrain(i,j)/=7) iso_max=max(iso_max,istream_order(i1-1,j1+1))
    if(idrain(i1-1,j1  )==8 .and. idrain(i,j)/=8) iso_max=max(iso_max,istream_order(i1-1,j1  ))
    if(idrain(i1-1,j1-1)==1 .and. idrain(i,j)/=1) iso_max=max(iso_max,istream_order(i1-1,j1-1))
    if(idrain(i1  ,j1-1)==2 .and. idrain(i,j)/=2) iso_max=max(iso_max,istream_order(i1  ,j1-1))
    if(idrain(i1+1,j1-1)==3 .and. idrain(i,j)/=3) iso_max=max(iso_max,istream_order(i1+1,j1-1))
    if(idrain(i1+1,j1  )==4 .and. idrain(i,j)/=4) iso_max=max(iso_max,istream_order(i1+1,j1  ))
    !if (iso_max>0) write(*,*)'iso_max',iso_max
    if(istream_order(i,j)<iso_max)then
      istream_order(i1,j1)=iso_max
    else if(istream_order(i,j)==iso_max)then
      istream_order(i1,j1)=iso_max+1
    else if(istream_order(i,j)>iso_max)then
      istream_order(i1,j1)=istream_order(i,j)
    else
      write(*,*)'from',istream_order(i,j)
      write(*,*)' to ',istream_order(i1,j1)
    end if
    i2=i
    j2=j
    i=i1
    j=j1
    !write(*,*)'iteration =',itr
   ! itr=itr+1
  end do
  !close(100)
end do
close(20,status='delete')
!
fname='streamOrder.asc'
call w_iasc(20,fname,istream_order(:,:),iend,jend,xllcorner,yllcorner,cellsize,0)
!fname='IO/iriver.asc'
!call w_iasc(20,fname,iriver(:,:),iend,jend,xllcorner,yllcorner,cellsize,0)
n_1d=count(istream_order(:,:)/=0)
!write(*,*)'n_1d=',n_1d
allocate(i_1d(n_1d),j_1d(n_1d),in_1d(n_1d),jn_1d(n_1d),ip_1d(n_1d),jp_1d(n_1d),idrain_1d(n_1d),iso_1d(n_1d),iacc_1d(n_1d))
allocate(x_1d(n_1d),y_1d(n_1d),z_1d(n_1d))!,acc_1d(n_1d))
i_1d(1:n_1d)=pack(i_2d(:,:),istream_order(:,:)/=0)
j_1d(1:n_1d)=pack(j_2d(:,:),istream_order(:,:)/=0)
in_1d(1:n_1d)=pack(i1_2d(:,:),istream_order(:,:)/=0)
jn_1d(1:n_1d)=pack(j1_2d(:,:),istream_order(:,:)/=0)
ip_1d(1:n_1d)=pack(i2_2d(:,:),istream_order(:,:)/=0)
jp_1d(1:n_1d)=pack(j2_2d(:,:),istream_order(:,:)/=0)
!
iacc_1d(1:n_1d)=pack(iacc(:,:),istream_order(:,:)/=0)
idrain_1d(1:n_1d)=pack(idrain(:,:),istream_order(:,:)/=0)
iso_1d(1:n_1d)=pack(istream_order(:,:),istream_order(:,:)/=0)
!
x_1d(1:n_1d)=pack(x_2d(:,:),istream_order(:,:)/=0)
y_1d(1:n_1d)=pack(y_2d(:,:),istream_order(:,:)/=0)
z_1d(1:n_1d)=pack(z_2d(:,:),istream_order(:,:)/=0)
!open(20,file='IO/river_struc.txt')
!write(20,'(a)')'x,y,z,i,j,i-1,j-1,i+1,j+1,km2,iacc,drain,stream_order'
!do i=1,n_1d
!  !write(*,'(2f15.3,6i5)')x_1d(i),y_1d(i),i_1d(i),j_1d(i),in_1d(i),jn_1d(i),ip_1d(i),jp_1d(i)
!  ba=cellsize*cellsize*float(iacc_1d(i))/1000000.d0
!  if(i_1d(i)==0 .or. j_1d(i)==0)then
!    write(*,'(a,10i5)')' i=0 or j=0 @',i_1d(i),j_1d(i),ip_1d(i),jp_1d(i),in_1d(i),jn_1d(i)
!  end if
!  write(20,'(3(f15.3,a),6(i5,a),(f10.5,a),(i10,a),2(i5,a))')x_1d(i),',',y_1d(i),',',z_1d(i),',',i_1d(i),',',j_1d(i),',',  &
!  ip_1d(i),',',jp_1d(i),',',in_1d(i),',',jn_1d(i),',',ba,',',iacc_1d(i),',',idrain_1d(i),',',iso_1d(i)
!end do
!close(20)
allocate(id_from(n_1d,3),id_to(n_1d),idx(3),dx(n_1d),b_1d(n_1d),d_1d(n_1d))!,acc_1d(n_1d))
! Width & depth Tanaka & Sayama 2018 DPRI
! Width
wm=4.73d0
wp=0.58d0
wmin=5.d0
! Depth
dm=1.57d0
dp=0.33d0
dmin=1.d0
!
mn=0.1d0 ! Manning coefficient
dsoil=1.d0
!
do i=1,n_1d
  id_from(i,:)=0
  itr=0
  do j=1,n_1d
    if(in_1d(j)==i_1d(i) .and. jn_1d(j)==j_1d(i))then
      itr=itr+1
      id_from(i,itr)=j
    end if
  end do
  idx(1:itr)=id_from(i,1:itr)
  iexit=0
  do j=1,itr
    mask=iso_1d(idx(j))==maxval(iso_1d(idx(1:itr)))
    if(mask .and. iacc_1d(idx(j))==maxval(iacc_1d(idx(1:itr))) .or. &
       mask .and. count(iso_1d(idx(1:itr))==maxval(iso_1d(idx(1:itr))))==1)then
      id_from(i,1:itr)=cshift(idx(1:itr),j-1)
      exit
    end if
  end do
end do
!
do i=1,n_1d
  idx(:)=0
  iexit=0
  do j=1,n_1d
    do k=1,3
      if(id_from(j,k)==i .and. minval(id_from(j,1:k))>=1)then
        idx(1)=j
        iexit=1
        exit
      end if
    end do
    if(iexit==1)exit
  end do
  !
  do j=1,n_1d
    if(ip_1d(j)==i_1d(i) .and. jp_1d(j)==j_1d(i))then
      idx(2)=j
    end if
  end do
  !
  id_to(i)=maxval(idx(:))
  if(id_to(i)==0)then
    dx(i)=-999
  else
    dx0=x_1d(id_to(i))-x_1d(i)
    dy0=y_1d(id_to(i))-y_1d(i)
    dx(i)=(dx0*dx0+dy0*dy0)**0.5d0
  end if
  ba=cellsize*cellsize*float(iacc_1d(i))/1000000.d0
  b_1d(i)=max(wmin,wm*ba**wp)
  d_1d(i)=max(dmin,dm*ba**dp)
end do
!
open(20,file='streamConfiguration_inRR.txt')
write(20,'(i5,a)')n_1d,',# of node'
!write(20,'(a)')'x,y,z,b,dep,dx,id,id_to,id_from1,id_from2,id_from3,iacc,i,j'
write(20,'(a)')'id,id_to,id_from1,id_from2,id_from3,z_2d,z_1d,zs_1d,dx,b,n,iacc,x_2d,y_2d,i_2d,j_2d'
do i=1,n_1d
  if(id_to(i)==0)then
    if(mod(idrain_1d(i)*(-45.d0)+90.d0,90.d0)==0)then
      !dx(i)=1.d0
      dx(i)=maxval(dx)/2.d0**0.5d0
    else
      !dx(i)=2.d0**0.5d0
      dx(i)=maxval(dx)
    end if
  end if
  write(20,'(10(i5,:'',''))',advance='no')i,id_to(i),id_from(i,1:3)
  write(20,'(a,6(f10.3,'',''))',advance='no')',',z_1d(i),z_1d(i)-d_1d(i),z_1d(i)-d_1d(i)-dsoil,dx(i),b_1d(i),mn
  !write(20,'((i10,:'',''),(f15.5,:'',''))',advance='no')ifld,acc(i)-sum(acc(id_from(i,1:3)))
  !write(20,'((i10,:'',''),(f15.5,:'',''))',advance='no')ifld,acc(i)
  write(20,'((i5,:'',''))',advance='no')int(iacc_1d(i))!-sum(acc(id_from(i,1:3)))
  write(20,'(a,2(f15.5,:'',''),2(i5,:'',''))')',',x_1d(i),y_1d(i),i_1d(i),j_1d(i)
end do
close(20)
!
open(20,file='streamConfiguration_Checker.txt')
write(20,'(a)')'x,y,b,i,i_to,dx,deg,grad'
do i=1,n_1d
  idx(1)=id_from(i,1)
  idx(2)=id_to(i)
  if(id_to(i)==0)then
    cycle
    deg=idrain_1d(i)*(-45.d0)+90.d0
    x_vec=x_1d(i)+dx(i)*cos(idrain_1d(i)*45.d0*deg2rad)*0.5d0
    y_vec=y_1d(i)+dx(i)*sin(idrain_1d(i)*45.d0*deg2rad)*0.5d0
    grad=0.d0
  else
    x_vec=(x_1d(i)+x_1d(idx(2)))*0.5d0
    y_vec=(y_1d(i)+y_1d(idx(2)))*0.5d0
    deg=idrain_1d(i)*(-45.d0)+90.d0
    dl=((x_1d(i)-x_1d(idx(2)))**2.0d0+(y_1d(i)-y_1d(idx(2)))**2.0d0)**0.5d0
    grad=(z_1d(idx(2))-z_1d(i))/dl
      !write(20,'(2f10.3,2i5,2f10.3)')x_vec,y_vec,i,idx(2),dx(i),deg
    write(20,'(3(f15.3,'',''))',advance='no')x_vec,y_vec,b_1d(i)
    write(20,'(2(i5,'',''))',advance='no')i,idx(2)
    write(20,'(3(f15.3,:'',''))')dx(i),deg,atan(-grad)*rad2deg
  end if
  !
  !write(20,'(2f10.3,2i5,2f10.3)')x_vec,y_vec,i,idx(2),dx(i),deg
  !write(20,'(3(f15.3,'',''))',advance='no')x_vec,y_vec,b_1d(i)
  !write(20,'(2(i5,'',''))',advance='no')i,idx(2)
  !write(20,'(3(f15.3,:'',''))')dx(i),deg,atan(-grad)*rad2deg
  !
end do
close(20)
!
open(20,file='dir_deg.txt')
write(20,*)'x,y,i,j,dx,deg'
do i=1,iend
  do j=1,jend
    if(ibasin(i,j)/=0)then
      if(mod(idrain(i,j),2)==0)then
        dl=cellsize
      else
        dl=cellsize*2.d0**0.5d0
      end if
      x=xllcorner+cellsize*(0.5+i-1)
      y=yllcorner+cellsize*(0.5+j-1)
      deg=idrain(i,j)*(45.d0)
      x=x+dl*0.5d0*cos(deg2rad*deg)
      y=y+dl*0.5d0*sin(deg2rad*deg)
      deg=idrain(i,j)*(-45.d0)+90.d0
      write(20,'(2(f15.3,'',''),2(i5,:'',''),2(f10.3,:'',''))')x,y,i,j,dl,deg
    end if
  end do
end do
close(20)
!
write(*,'(a)')'--- nomal end ---'
stop
end
!
!
!
subroutine f_next(iend,jend,idrain,i,j,i1,j1)
implicit none
integer :: i,j,iend,jend,i1,j1
integer :: idrain(iend,jend)
!
! nw n ne ! 3  2  1
!  \ | /  !  \ | /
! w -+- e ! 4 -+- 8
!  / | \　　!  / | \
! sw s se　! 5  6 　7　　　　
!idrain(i,j)=abs(idrain(i,j))
if(     idrain(i,j)==1)then
  i1=i+1
  j1=j+1
else if(idrain(i,j)==2)then
  i1=i
  j1=j+1
else if(idrain(i,j)==3)then
  i1=i-1
  j1=j+1
else if(idrain(i,j)==4)then
  i1=i-1
  j1=j
else if(idrain(i,j)==5)then
  i1=i-1
  j1=j-1
else if(idrain(i,j)==6)then
  i1=i
  j1=j-1
else if(idrain(i,j)==7)then
  i1=i+1
  j1=j-1
else if(idrain(i,j)==8)then
  i1=i+1
  j1=j
!else if(idrain(i,j)<=-1 .or. iriver(i,j)==0)then
!  write(*,'(a,4i5)')'river end @ i,j,dir,iriv=',i,j,idrain(i,j),iriver(i,j)
!  itmp=1
else if(i1<1 .or. i1>iend .or. j1<1 .or. j1>jend)then
!  write(*,'(a,4i5)')'river end @ i,j,dir,iriv=',i,j,idrain(i,j)
  i1=-1
  j1=-1
  !itmp=1
else
  write(*,'(a,4i5)')'river course? @ i,j,dir,iriv=',i,j,idrain(i,j) !,iriver(i,j)
  write(*,*)'Don''t come here?'
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
