!
! Written by "Ryusuke NUMATA" <ryusuke.numata@gmail.com>
!                              https://rnumata.org
!         in 28 May 2016.
!
! $Id$
!
! see https://rnumata.org/simul/mpiex1
program diffusion_mpi

  use mpi

  implicit none

  integer :: NX0=512               ! grid points [input]
  integer :: NT=100000             ! iterations [input]
  integer :: NOUT=1000             ! output size [input]
  integer, parameter :: NBX=1      ! boundary grid points
  real, parameter :: DC=1.         ! diffusion coefficient
  real, parameter :: LENGTH=2.     ! domain size [-LENGTH/2:LENGTH/2]
  integer :: nx
  integer :: ix,it
  real, allocatable :: f(:),df(:),x(:)
  real, allocatable :: fp(:),xp(:)
  real, allocatable :: fall(:),xall(:)
  real :: dx
  real :: t,dt
  real, parameter :: pi=4.*atan(1.)
  real :: total,total0

  include 'diffuse.h'

  integer, parameter :: rootproc=0
  integer :: myid=rootproc,nps=1,ip,nrange
  integer, allocatable :: npmin(:),npmax(:),npbmin(:),npbmax(:)
  integer :: ierr
  integer :: p_send,p_recv
  integer :: mpirealtype,tag=0,status(mpi_status_size)
  integer, allocatable :: gv_count(:), gv_dspl(:)

  namelist /input/ NX0,NT,NOUT

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,myid,ierr)
  call mpi_comm_size(mpi_comm_world,nps,ierr)

  mpirealtype=mpi_real

  allocate(npmin(0:nps-1),npmax(0:nps-1))
  allocate(npbmin(0:nps-1),npbmax(0:nps-1))
  allocate(gv_count(0:nps-1),gv_dspl(0:nps-1))
  
  ! input
  if(myid==rootproc) read(5,input)
  call mpi_bcast(NX0, 1,mpi_integer,rootproc,mpi_comm_world,ierr)
  call mpi_bcast(NT,  1,mpi_integer,rootproc,mpi_comm_world,ierr)
  call mpi_bcast(NOUT,1,mpi_integer,rootproc,mpi_comm_world,ierr)

  nrange=NX0/nps
  do ip=0,nps-1
     npmin(ip)=nrange*ip+1  +NBX
     npmax(ip)=nrange*(ip+1)+NBX
  end do
  npmax(nps-1)=NX0+NBX
  npbmin(0:nps-1)=npmin(0:nps-1)-NBX
  npbmax(0:nps-1)=npmax(0:nps-1)+NBX

  ! setup gatherv
  if (nps==1) then
     gv_count(0)=npbmax(0)-npbmin(0)+1
     gv_dspl(0)=npbmin(0)-1
  else
     gv_count(0)=npmax(0)-npbmin(0)+1
     gv_dspl(0)=npbmin(0)-1
     do ip=1,nps-2
        gv_count(ip)=npmax(ip)-npmin(ip)+1
        gv_dspl(ip)=npmin(ip)-1
     end do
     gv_count(nps-1)=npbmax(nps-1)-npmin(nps-1)+1
     gv_dspl(nps-1)=npmin(nps-1)-1
  end if
  
  nx=NX0+2*NBX

  allocate(f(npbmin(myid):npbmax(myid)),df(npbmin(myid):npbmax(myid)),x(npbmin(myid):npbmax(myid)))
  if (myid==rootproc) then
     allocate(fp(0:NX0),xp(0:NX0))
     allocate(fall(nx),xall(nx))
  end if

  dx = LENGTH/NX0

  do ix = npbmin(myid), npbmax(myid)
     x(ix) = dx*(ix-NBX-.5) - .5*LENGTH
  end do

  call mpi_gatherv( &
       & x(gv_dspl(myid)+1),gv_count(myid),mpirealtype, &
       & xall,gv_count,gv_dspl,mpirealtype, &
       & rootproc,mpi_comm_world,ierr)
  
  dt = .5 * dx * dx / DC
  t = 0.

  !
  ! initialize array
  !
  f(npbmin(myid):npbmax(myid))=0.
  df(npbmin(myid):npbmax(myid))=0.

  !
  ! initial condition
  !
  f(npbmin(myid):npbmax(myid)) = .5*(cos(2.*pi*x(npbmin(myid):npbmax(myid))/LENGTH)+1.)

  call mpi_gatherv( &
       & f(gv_dspl(myid)+1),gv_count(myid),mpirealtype, &
       & fall,gv_count,gv_dspl,mpirealtype, &
       & rootproc,mpi_comm_world,ierr)
  
  if(myid==rootproc) then
     !
     ! write parameter for gnuplot
     !
     open(20,file=GPFILE)
     write(20,*) "nx0 = ",nx0
     write(20,*) "nt = ",nt
     write(20,*) "nout = ",nout
     write(20,*) "dt = ",dt
     close(20)
  
     open(10,file=DATAFILE)

     !
     ! write annotations
     !
1000 format( &
          "# Parameters" / &
          "#   NX0 = ",i0,", NT = ",i0,", NBX = ",i0,", NOUT = ",i0, / &
          "#   DC = ",f10.5," LENGTH = ",f10.5 / &
          "#" )
     write(10,1000) NX0, NT, NBX, NOUT, DC, LENGTH
     write(10,'("#",a24,2(1x,a24))')'time','x','f(x)'

     !
     ! write initial data
     !
     xp(0:NX0)=.5*(xall(NBX:NX0+NBX)+xall(NBX+1:NX0+NBX+1))
     fp(0:NX0)=.5*(fall(NBX:NX0+NBX)+fall(NBX+1:NX0+NBX+1))
     do ix=0,NX0
        write(10,'(3(1x,e24.12))') t,xp(ix),fp(ix)
     end do
     write(10,*)
     
     total=sum(fall(1+NBX:NX0+NBX))*dx
     total0=total
     write(6,'(3(a,f12.4),a)') &
          & 'time = ',t, &
          & ', total = ',total, &
          & ', difference = ',(total0-total)/total0*100,' [%]'
  end if

  !
  ! main routine
  !
  do it = 1, NT

     !
     ! advance
     !
     df(npmin(myid):npmax(myid)) = DC*dt/dx**2 * &
          & (    f(npmin(myid)+1:npmax(myid)+1) &
          & - 2.*f(npmin(myid)  :npmax(myid)) &
          & +    f(npmin(myid)-1:npmax(myid)-1) )
     f(npmin(myid):npmax(myid)) = &
          & f(npmin(myid):npmax(myid)) + df(npmin(myid):npmax(myid))

     !
     ! boundary condition
     !
     if(nps==1) then

        f(npbmin(myid):npmin(myid)-1)=f(npmax(myid)-NBX+1:npmax(myid))
        f(npmax(myid)+1:npbmax(myid))=f(npmin(myid):npmin(myid)+NBX-1)

     else

        !
        ! Periodic boundary
        ! Processor 0 sends f(NBX+1:NBX+NBX) and 
        ! Processor nps-1 receives/stores it in f(NX0+NBX+1:NX0+NBX+NBX)
        !
        p_send=mpi_proc_null
        p_recv=mpi_proc_null
        if(myid==0)     p_send=nps-1
        if(myid==nps-1) p_recv=0
        call mpi_send(f(1+NBX),NBX,mpirealtype, &
             & p_send,tag,mpi_comm_world,ierr)
        call mpi_recv(f(NX0+NBX+1),NBX,mpirealtype, &
             & p_recv,tag,mpi_comm_world,status,ierr)

        !
        ! Periodic boundary
        ! Processor nps-1 sends f(NX0+1:NX0+NBX) and
        ! Processor 0 receives/stores it in f(1:NBX)
        !
        p_send=mpi_proc_null
        p_recv=mpi_proc_null
        if(myid==0)     p_recv=nps-1
        if(myid==nps-1) p_send=0
        call mpi_send(f(NX0+1),NBX,mpirealtype, &
             & p_send,tag,mpi_comm_world,ierr)
        call mpi_recv(f(1),NBX,mpirealtype, &
             & p_recv,tag,mpi_comm_world,status,ierr)

        !
        ! Shift communications
        ! send/receive boundary data to right/from left
        !
        p_send=myid+1
        p_recv=myid-1
        if(myid==0)     p_recv=mpi_proc_null
        if(myid==nps-1) p_send=mpi_proc_null
        call mpi_send(f(npmax(myid)-NBX+1),NBX,mpirealtype, &
             & p_send,tag,mpi_comm_world,ierr)
        call mpi_recv(f(npmin(myid)-NBX),NBX,mpirealtype, &
             & p_recv,tag,mpi_comm_world,status,ierr)

        !
        ! Shift communications
        ! send/receive boundary data to left/from right
        !
        p_send=myid-1
        p_recv=myid+1
        if(myid==0)     p_send=mpi_proc_null
        if(myid==nps-1) p_recv=mpi_proc_null
        call mpi_send(f(npmin(myid)),NBX,mpirealtype, &
             & p_send,tag,mpi_comm_world,ierr)
        call mpi_recv(f(npmax(myid)+1),NBX,mpirealtype, &
             & p_recv,tag,mpi_comm_world,status,ierr)

     end if

     call mpi_gatherv( &
          & f(gv_dspl(myid)+1),gv_count(myid),mpirealtype, &
          & fall,gv_count,gv_dspl,mpirealtype, &
          & rootproc,mpi_comm_world,ierr)
     
     !
     ! output data
     !
     t = t + dt
     if(mod(it,NOUT)==0 .and. myid==rootproc) then
        fp(0:NX0)= .5*(fall(NBX:NX0+NBX)+fall(NBX+1:NX0+NBX+1))

        do ix=0,NX0
           write(10,'(3(1x,e24.12))') t,xp(ix),fp(ix)
        end do
        write(10,*)

        total=sum(fall(1+NBX:NX0+NBX))*dx
        write(6,'(3(a,f12.4),a)') &
             & 'time = ',t, &
             & ', total = ',total, &
             & ', difference = ',(total0-total)/total0*100.,' [%]'
     end if

  end do
  close(10)

  call mpi_finalize(ierr)

  stop

end program diffusion_mpi
