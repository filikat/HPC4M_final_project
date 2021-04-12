program final
implicit none
include 'mpif.h'
integer,parameter:: n=400
real*8:: To,Ti,T(0:n,0:n),res,Told,t_start,t_end,tolerance
real*8:: t_int1,t_int2
integer:: M(0:n,0:n),i,j,k,iter,holetype
integer:: comm,rank,ierr,nproc
integer:: total_nodes,column_nodes(1:n-1),nodes4proc,col
integer,allocatable::col4proc(:),Mproc(:,:)
real*8,allocatable::Tproc(:,:),resother(:)
integer:: start,sumcol,finish,dimsend,maxiter
logical:: cont

comm=MPI_COMM_WORLD
call mpi_init(ierr)
call mpi_comm_rank(comm,rank,ierr)
call mpi_comm_size(comm,nproc,ierr)


call mpi_barrier(comm,ierr)
t_start=mpi_wtime()

!inner/outer temperature
To=0.d0
Ti=1.d0

!parameters for the iteration
maxiter=50000
tolerance=1e-8

holetype=3
!holetype=  1-rectangular  2-circular
!           3-odd shape    4-concave


!Create information about the domain in node 0
!and send it to the other processors
if (rank==0) then

  allocate(col4proc(1:nproc),resother(1:nproc))

  ! Create domain:
  ! %% information matrix
  M=1
  ! % border
  M(0,:)=2
  M(n,:)=2
  M(:,0)=2
  M(:,n)=2
  ! %% hole
  call hole(M,n,holetype)

  total_nodes=sum(M(1:n-1,1:n-1))           !total nodes that require computation
  column_nodes=sum(M(1:n-1,1:n-1),dim=2)    !nodes per column
  nodes4proc=int(total_nodes/nproc)         !~ nodes per processor


  !sum nodes in columns until they reach 90% of nodes4proc
  start=1
  do i=1,nproc-1
    sumcol=0
    j=0
    do while (sumcol<9*nodes4proc/10)
      sumcol=sumcol+column_nodes(start+j)
      j=j+1
    enddo
    start=start+j
    col4proc(i)=j
  enddo
  col4proc(nproc)=n-1-sum(col4proc(1:nproc-1))

  !split information M into pieces and send them to the respective processor
  allocate(Mproc(0:n,0:col4proc(1)+1))
  Mproc=M(0:n,0:col4proc(1)+1)
  start=col4proc(1)
  do i=2,nproc
    dimsend=(n+1)*(col4proc(i)+2)
    finish=start+col4proc(i)+1
    call mpi_send(M(0:n,start:finish),dimsend,mpi_int,i-1,i-1,comm,ierr)
    start=start+col4proc(i)
  enddo

endif

call mpi_barrier(comm,ierr)
t_int1=mpi_wtime()

if (rank.ne.0) then
  allocate(col4proc(1:nproc))
endif

!nodes receive information about how many columns they have been assigned
call mpi_bcast(col4proc,nproc,mpi_int,0,comm,ierr)

!nodes receive their share of the domain
if (rank.ne.0) then
  allocate(Mproc(0:n,0:col4proc(rank+1)+1))
  call mpi_recv(Mproc,(n+1)*(col4proc(rank+1)+2),mpi_int,0,rank,comm,mpi_status_ignore,ierr)
endif

!Now all nodes have the respective matrix Mproc, from which
!they create the matrix Tproc with the temperature data
col=col4proc(rank+1)
allocate(Tproc(0:n,0:col+1))
! %% temperature matrix
Tproc=(Ti+To)/2
where (Mproc==2)
  Tproc=To
elsewhere (Mproc==0)
  Tproc=Ti
end where

! %% iteration
res=1.d0
iter=0
cont=.true.

!iteration continues until it reaches maximum iterations or
!until processor 0 gives the order to stop (using cont)
do while ((iter<maxiter).and.(cont))
  iter=iter+1
  if (res>tolerance) then
    res=0.d0
    do i=1,n-1
      do j=1,col
        if (Mproc(i,j)==1) then      !perform computation only for nodes that are internal
          Told=Tproc(i,j)
          Tproc(i,j)=(Tproc(i,j-1)+Tproc(i,j+1)+Tproc(i-1,j)+Tproc(i+1,j))/4
          res=max(res,abs(Told-Tproc(i,j)))     !res is the maximum residual among the nodes in the domain
        endif
      enddo
    enddo
  endif

  !perform halo swapping
  if (rank==0) then
    call mpi_send(Tproc(0:n,col),n+1,mpi_real8,rank+1,1,comm,ierr)
    call mpi_recv(Tproc(0:n,col+1),n+1,mpi_real8,rank+1,1,comm,mpi_status_ignore,ierr)
  elseif (rank==nproc-1) then
    call mpi_recv(Tproc(0:n,0),n+1,mpi_real8,rank-1,1,comm,mpi_status_ignore,ierr)
    call mpi_send(Tproc(0:n,1),n+1,mpi_real8,rank-1,1,comm,ierr)
  else
    call mpi_recv(Tproc(0:n,0),n+1,mpi_real8,rank-1,1,comm,mpi_status_ignore,ierr)
    call mpi_send(Tproc(0:n,1),n+1,mpi_real8,rank-1,1,comm,ierr)
    call mpi_send(Tproc(0:n,col),n+1,mpi_real8,rank+1,1,comm,ierr)
    call mpi_recv(Tproc(0:n,col+1),n+1,mpi_real8,rank+1,1,comm,mpi_status_ignore,ierr)
  endif

  !processor 0 collects the residuals from every processor and
  !decides if the computation needs to be stopped
  call mpi_gather(res,1,mpi_real8,resother,1,mpi_real8,0,comm,ierr)
  if ((rank==0).and.(maxval(resother)<tolerance)) then
      cont=.false.
      deallocate(resother)
  endif
  call mpi_bcast(cont,1,mpi_logical,0,comm,ierr)

enddo

print *, 'node ',rank,', final residual: ',res

call mpi_barrier(comm,ierr)
t_int2=mpi_wtime()

!Iteration finished, processor 0 now builds a single matrix
!with the information in each processor

!each processor sends the solution to proc 0.
!they send only the internal columns, not the first and last one.
if (rank.ne.0) then
  call mpi_send(Tproc(0:n,1:col),(n+1)*col,mpi_real8,0,rank,comm,ierr)
endif

!processor 0 puts the bits of solution together and writes to file
if (rank==0) then
  T(0:n,1:col)=Tproc(0:n,1:col)
  start=col+1
  do i=1,nproc-1
    finish=start-1+col4proc(i+1)
    call mpi_recv(T(0:n,start:finish),(n+1)*col4proc(i+1),mpi_real8,i,i,comm,mpi_status_ignore,ierr)
    start=start+col4proc(i+1)
  enddo

  open(10,file='data.dat')
  do i=0,n
    write(10,*) T(i,:)
  enddo

endif

call mpi_barrier(comm,ierr)
t_end=mpi_wtime()
if (rank==0) then
  print '(A,I15)',   ' Hole type      : ', holetype
  print '(A,I15)',   ' Subdivisions   : ', n
  print '(A,I15)',   ' Iterations     : ', iter
  print '(A,ES15.2)',' Final residual : ', res
  print '(A,I15)',   ' Processors     : ',nproc
  print '(A,F15.4)', ' Total Time     : ',t_end-t_start
  print '(A,F15.4)', ' Preprocessing  : ',t_int1-t_start
  print '(A,F15.4)', ' Postprocessing : ',t_end-t_int2
endif

deallocate(col4proc,Mproc,Tproc)
call mpi_finalize(ierr)

contains


  subroutine hole(M,n,holetype)
    implicit none
    integer,intent(in):: n,holetype
    integer,intent(inout)::M(0:n,0:n)
    integer:: i,j

  if (holetype==1) then
    !rectangular hole
    M(int(n/2-n/8):int(n/2+n/3),int(n/2-n/12):int(n/2+n/5))=0

  elseif (holetype==2) then
    !circular hole
    do i=0,n
      do j=0,n
        if ((i-n*9/20)**2+(j-n*13/20)**2<(n/5)**2) then
          M(i,j)=0
        endif
      enddo
    enddo

  elseif (holetype==3) then
    !strange hole
    M(n/4:3*n/8,3*n/16:5*n/16)=0
    M(n/4:5*n/8,5*n/16:3*n/8)=0
    M(3*n/8:3*n/4,3*n/8:5*n/8)=0
    M(5*n/8:3*n/4,5*n/8:7*n/8)=0

  elseif (holetype==4) then
    !concave hole
    M(3*n/16:3*n/8,9*n/16:9*n/16+10)=0
    M(3*n/8:3*n/8+10,3*n/8:11*n/16)=0
    M(3*n/8+10:5*n/8,3*n/8:3*n/8+10)=0
    M(3*n/8+10:9*n/16,11*n/16-10:11*n/16)=0


  endif

  end subroutine hole
end program final
