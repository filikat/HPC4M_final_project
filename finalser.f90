program final
implicit none
integer,parameter:: n=100
real*8::dx,To,Ti,T(0:n,0:n),res,Told
integer:: M(0:n,0:n),i,j,k,iter,type

type=4
!type=  1-rectangular  2-circular
!       3-odd shape    4-concave

To=0.d0
Ti=1.d0
dx=1.d0/n

! Create domain:
! %% information matrix
M=1
! % border
M(0,:)=0
M(n,:)=0
M(:,0)=0
M(:,n)=0
! %% hole
call hole(M,n,type)

! %% temperature matrix
T=(Ti+To)/2
where (M==0)
  T=To
elsewhere (M==2)
  T=Ti
end where

! %% iteration
res=1.d0
iter=0
do while ((res>1e-6).and.(iter<20000))
  res=0.d0
  iter=iter+1
  do i=1,n-1
    do j=1,n-1
      if (M(i,j)==1) then
        Told=T(i,j)
        T(i,j)=(T(i,j-1)+T(i,j+1)+T(i-1,j)+T(i+1,j))/4
        res=max(res,abs(Told-T(i,j)))
      endif
    enddo
  enddo
enddo

print '(A,I15)', ' Hole type      : ', type
print '(A,I15)',' Subdivisions   : ', n
print '(A,I15)',' Iterations     : ', iter
print '(A,ES15.2)',' Final residual : ', res

open(10,file='data.dat')
do i=0,n
  write(10,*) T(i,:)
enddo


contains


  subroutine hole(M,n,type)
    implicit none
    integer,intent(in):: n,type
    integer,intent(inout)::M(0:n,0:n)
    integer:: i,j

  if (type==1) then
    !rectangular hole
    !M(n/2-15:n/2+3,n/2-3:n/2+3)=2
    M(n/2-15:n/2+3,n/2-3:n/2+3)=2

  elseif (type==2) then
    !circular hole
    do i=0,n
      do j=0,n
        if ((i-n/2)**2+(j-n/2)**2<(n/10)**2) then
          M(i,j)=2
        endif
      enddo
    enddo

  elseif (type==3) then
    !strange hole
    M(60:70,20:40)=2
    M(50:60,20:50)=2
    M(40:50,30:60)=2
    M(30:40,40:70)=2
    M(20:30,50:70)=2

  elseif (type==4) then
    !concave hole
    M(20:30,30:60)=2
    M(30:60,30:40)=2
    M(60:70,30:60)=2
    M(40:50,20:30)=2
  endif

  end subroutine hole
end program final
