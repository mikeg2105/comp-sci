
!
! function routines for laxsolver
!

module laxsolver


implicit none

public

public :: halfstep, flux, sadvance


! momentum, density, energy, pressure - half-step values
double precision ,allocatable :: mxhm(:),rhohm(:),ehm(:),phm(:)

! momentum, density, energy, pressure - half-step values
double precision ,allocatable :: mxhp(:),rhohp(:),ehp(:),php(:)

! momentum, density, energy, pressure - half-step values
double precision,allocatable :: fmx(:),frho(:),fe(:),fp(:)



! constants to use

double precision :: mmx, mrho, me, mp







contains

subroutine initsolverarrays(ni,nt)
    integer, intent(in) :: ni, nt


    allocate (fmx(1:ni))
    allocate (frho(1:ni))
    allocate (fe(1:ni))
    allocate (fp(1:ni))

    allocate (mxhm(1:ni))
    allocate (rhohm(1:ni))
    allocate (ehm(1:ni))
    allocate (phm(1:ni))

    allocate (mxhp(1:ni))
    allocate (rhohp(1:ni))
    allocate (ehp(1:ni))
    allocate (php(1:ni))




    mxhm=0
    rhohm=0
    ehm=0
    phm=0

    mxhp=0
    rhohp=0
    ehp=0
    php=0

    fmx=0
    frho=0
    fe=0
    fp=0

endsubroutine





!
!    subroutine halfstep(timestep, deltax, denscurrent, dens, densh)
    subroutine halfstep(ni, nt, it, deltat, deltax, flux, field, fieldhdp,fieldhdm)
            integer, intent(in) :: it, ni, nt
            double precision, intent(in) :: deltat, deltax
       		double precision, dimension(ni), intent(in) :: flux
       		double precision, dimension(ni,nt), intent(in) :: field
       		double precision, dimension(ni), intent(inout)  :: fieldhdp,fieldhdm
            double precision :: deltav
            integer :: i

            deltav=deltat/deltax

            do i=1,ni-1
                fieldhdp(i)=-(deltav/2.0d0)*(flux(i+1)-flux(i))+((field(i,it)+field(i+1,it))/2.0d0)
            end do

            do i=2,ni
                fieldhdm(i)=-(deltav/2.0d0)*(flux(i)-flux(i-1))+((field(i-1,it)+field(i,it))/2.0d0)
            end do



    endsubroutine


    subroutine flux(ifield, ni,it, mx, rho, e, p, fflux, fgamma )
            integer, intent(in) :: ifield, it, ni
       		double precision, dimension(ni), intent(inout) :: fflux
       		double precision, dimension(ni), intent(in) :: mx, rho, e, p
            double precision :: fgamma
            integer :: i

!mx:1, rho:2, e:3
!mx
        if(ifield.eq.1) then
            fflux=mx(:)
!rho
        elseif (ifield.eq.2) then
            fflux=-(mx(:)*mx(:)/rho(:))-(fgamma-1)*(rho(:)*e(:))
!e
        elseif (ifield.eq.3) then
            fflux=-mx(:)*e(:)/rho(:)
        else
            fflux=0
        endif



    endsubroutine

    !call sadvance(ifield, ni, nt, it, dt, dx, frho, rho, rhohm, rhohp, mxhm, mxhp, em, ep, fgamma)
    subroutine sadvance(ifield, ni, nt, it, deltat, deltax, field,  rhohm, rhohp, mxhm, mxhp, em, ep, fgamma)
            integer, intent(in) :: it, ni, nt, ifield
            double precision, intent(in) :: deltat, deltax, fgamma
       		double precision, dimension(ni,nt), intent(inout) :: field
       		double precision, dimension(ni), intent(in)  :: rhohm, rhohp, mxhm, mxhp, em, ep
            double precision :: deltav
            double precision, allocatable :: ghp(:),ghm(:)
            integer :: i

            allocate (ghp(1:ni))
            allocate (ghm(1:ni))

            deltav=deltat/deltax


            !calculate ghp and ghm for each field
            !mx
            if(ifield.eq.1) then
                do i=1,ni
                    ghp(i)=-(   (mxhp(i)*mxhp(i)/rhohp(i))+(fgamma-1)*ep(i)*rhohp(i))
                    ghm(i)=-(   (mxhm(i)*mxhm(i)/rhohm(i))+(fgamma-1)*em(i)*rhohm(i))
                enddo
            !rho
            elseif (ifield.eq.2) then
                    ghp=-mxhp
                    ghm=-mxhm
            !e
            elseif (ifield.eq.3) then
                do i=1,ni
                    ghp(i)=-((mxhp(i)*ep(i)/rhohp(i)))
                    ghm(i)=-((mxhm(i)*em(i)/rhohm(i)))
                enddo
            else
                    ghp=0
                    ghm=0
            endif

            field(1:ni,it+1)=field(1:ni,it)+deltav*(ghp-ghm)

    endsubroutine


end module laxsolver
