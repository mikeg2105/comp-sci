
!
! module handle input to and from files
!

module lparamdata


implicit none

public
character(len=20) :: filename

integer ni, nt, it, i
double precision xmax, dx
double precision dt, tmax, courant

! momentum, density, energy, pressure
double precision,allocatable :: mx(:,:),rho(:,:),e(:,:),p(:,:)


! time step and spatial interval
double precision,allocatable :: t(:),x(:)


! constants to use
double precision mx0, rho0, e0, p0, wavespeed, fgamma

!private

contains

subroutine initpatarrays(ni,nt)
    integer, intent(in) :: ni, nt

    allocate (mx(1:ni,1:nt))
    allocate (rho(1:ni,1:nt))
    allocate (e(1:ni,1:nt))
    allocate (p(1:ni,1:nt))

    allocate (t(1:nt))
    allocate (x(1:ni))


    ! Initialise configuration
    mx=0
    rho=0
    e=0
    p=0

    t=0
    x=0


endsubroutine



end module lparamdata
