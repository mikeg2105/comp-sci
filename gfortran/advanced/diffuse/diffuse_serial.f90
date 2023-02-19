!
! Written by "Ryusuke NUMATA" <ryusuke.numata@gmail.com>
!                              https://rnumata.org
!         in 28 May 2016.
!
! $Id$
!

program diffusion_serial

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
  real :: dx
  real :: t,dt
  real, parameter :: pi=4.*atan(1.)
  real :: total,total0

  include 'diffuse.h'

  namelist /input/ NX0,NT,NOUT

  ! input
  read(5,input)

  nx=NX0+2*NBX

  allocate(f(nx),df(nx),x(nx))
  allocate(fp(0:NX0),xp(0:NX0))

  dx = LENGTH/NX0

  do ix = 1, nx
     x(ix) = dx*(ix-NBX-.5) - .5*LENGTH
  end do
  xp(0:NX0)= .5*(x(NBX:NX0+NBX)+x(NBX+1:NX0+NBX+1))

  dt = .5 * dx * dx / DC
  t = 0.

  !
  ! initialize array
  !
  f(1:nx)=0.
  df(1:nx)=0.

  !
  ! initial condition
  !
  f(1:nx) = .5*(cos(2.*pi*x(1:nx)/LENGTH)+1.)
  fp(0:NX0)= .5*(f(NBX:NX0+NBX)+f(NBX+1:NX0+NBX+1))

  !
  ! write parameter for gnuplot
  !
  open(20,file=GPFILE)
  write(20,*) "nx0 = ",nx0
  write(20,*) "nt = ",nt
  write(20,*) "nout = ",nout
  write(20,*) "dt = ",dt
  close(20)
  
  !
  ! write annotations
  !
  open(10,file=DATAFILE)
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
  do ix=0,NX0
     write(10,'(3(1x,e24.12))') t,xp(ix),fp(ix)
  end do
  write(10,*)
  
  total=sum(f(1+NBX:NX0+NBX))*dx
  total0=total
  write(6,'(3(a,f12.4),a)') &
       & 'time = ',t, &
       & ', total = ',total, &
       & ', difference = ',(total0-total)/total0*100,' [%]'

  !
  ! main routine
  !
  do it = 1, NT

     !
     ! advance
     !
     df(1+NBX:NX0+NBX) = DC*dt/dx**2 * &
          & ( f(2+NBX:NX0+NBX+1) - 2. * f(1+NBX:NX0+NBX) + f(NBX:NX0+NBX-1) )
     f(1+NBX:NX0+NBX) = f(1+NBX:NX0+NBX) + df(1+NBX:NX0+NBX)

     !
     ! boundary condition
     !
     f(        1:        NBX)=f(NX0+    1:NX0+    NBX)
     f(NX0+NBX+1:NX0+NBX+NBX)=f(    NBX+1:    NBX+NBX)

     !
     ! output data
     !
     t = t + dt
     if(mod(it,NOUT)==0) then
        fp(0:NX0)= .5*(f(NBX:NX0+NBX)+f(NBX+1:NX0+NBX+1))
        do ix=0,NX0
           write(10,'(3(1x,e24.12))') t,xp(ix),fp(ix)
        end do
        write(10,*)

        total=sum(f(1+NBX:NX0+NBX))*dx
        write(6,'(3(a,f12.4),a)') &
             & 'time = ',t, &
             & ', total = ',total, &
             & ', difference = ',(total0-total)/total0*100.,' [%]'
     end if

  end do
  close(10)

  stop

end program diffusion_serial
