module statbalancemod
    implicit none
!mu=0.6d0;
!R=8.31e3;
! parrVALMc=rhoarrVALMc*TarrVALMc*R/mu

!rho, mom1, mom2, mom3, energy, b1, b2, b3,energyb,rhob,b1b,b2b,b3b
! p[*,*,*]=w[*,*,*,4]+w[*,*,*,8]
! p[*,*,*]=p[*,*,*]-(w[*,*,*,1]^2.0+w[*,*,*,2]^2.0+w[*,*,*,3]^2.0)/ $
!          (w[*,*,*,0]+w[*,*,*,9])/2.0
! p[*,*,*]=p[*,*,*]-((w[*,*,*,5]+w[*,*,*,10])^2.0+(w[*,*,*,6]+w[*,*,*,11])^2.0 $
!           +(w[*,*,*,7]+w[*,*,*,12])^2.0)/2.0
! p[*,*,*]=(gamma-1.d0)*p[*,*,*]


!compute correct pressure for gravitationally stratified atmosphere

!compute initial energy (at photosphere or temperature minimum)
!mu_thermal=0.6d0;
!R=8.31e3;

! temp*R*density/((mu_thermal))
!parrVALMc=rhoarrVALMc*TarrVALMc*R/mu
!iniene=6840.d0*8.31e3*(2.3409724e-09)/0.6d0/(eqpar(gamma_)-1.0)

! !iniene=731191.34d0*8.31e3*(1.1806882e-11)/0.6d0/(eqpar(gamma_)-1.0)
!
! !iniene=731191.34d0*8.31e3*(1.1790001e-11)/0.6d0/(eqpar(gamma_)-1.0)
!
! ! 1.6Mm
!
! !iniene=6840.d0*8.31e3*(2.3409724e-09)/0.6d0/(eqpar(gamma_)-1.0)
!iniene=6840.d0*R*(2.3409724e-09)/mu/(consts.fgamma-1.0);

!%
!% !iniene=6840.d0*8.31e3*(2.2139002e-09)/0.6d0/(eqpar(gamma_)-1.0)
!%
!% !iniene=731191.34d0*8.31e3*(4.5335481e-12)/0.6d0/(eqpar(gamma_)-1.0)



    real, parameter :: pi=4.*atan(1.0)
    real, parameter :: mu_mass=0.6d0
    real, parameter :: R=8.31e3

! adiabatic gas parameter and the magnetic permeablity
    real, parameter :: fgamma=1.66666667e0,mumag=4*pi/1.0e7
!   density kg/m^3
    real, parameter :: rho0=2.34d-4, p0=9228.6447

!   solar gravity m/s^2
    real, parameter :: gs=-274.0

!   temperature profile parameters (K))
    real, parameter :: Tch=8000, Tc=1.8d6

!   position of the transition zone (ytr) and width of transition zone (wtr) in metres
    real, parameter :: ytr=2.0d6, wtr=0.02d6


    public :: pi
    public :: mu_mass, R,fgamma,mumag
    public :: rho0, p0, gs
    public :: Tch, Tc, ytr, wtr
    public :: writefile, temp, hydropres, dens


private


contains

!compute temp at height using tanh function
real function temp( height )
    real, intent(in) :: height
    real :: tmptemp

    tmptemp=1+tanh((height-ytr)/wtr)

    temp=Tch+((Tc-Tch)/2.0)*tmptemp
end function

!compute pres
real function hydropres(heights, hindex, npoints, deltah)
    real, intent(in) :: deltah, heights(4096)
    integer, intent(in) :: hindex, npoints
    real :: psum,Hscale
    integer :: i

!    tmp0=temp(heights(1))
    psum=0.0

    if (hindex.eq.npoints) then
        Hscale=R*temp(heights(hindex))/(mu_mass*gs)
!        tmptemp=p0*tmp0/temp(heights(hindex))
        psum=psum+p0*exp(deltah/Hscale)
    elseif (hindex.lt.npoints) then
 !       tmptemp=p0*tmp0/temp(heights(hindex))
        do i=npoints,hindex,-1
            Hscale=R*temp(heights(i))/(mu_mass*gs)
            psum=psum+deltah/Hscale
        end do
        psum=p0*exp(psum)
    endif
    print*,'psum ',hindex,' ',psum
    hydropres=psum

end function

real function dens( height, spres )
    real, intent(in) :: height, spres
    real :: tmpresult

! parrVALMc=rhoarrVALMc*TarrVALMc*R/mu
    tmpresult= mu_mass*spres/(R*temp(height))
    dens=tmpresult

end function

subroutine writefile(height, dens, press, temp, nitems)
    implicit none
    integer, intent(in) :: nitems
    real, intent(in) :: height(nitems), dens(nitems), press(nitems), temp(nitems)

    integer :: i

    open(unit=9, file='atmos.txt')

    do i=1,nitems
        write(9,*) height(i),temp(i),dens(i), press(i)
200     format(F16.6,2X,F16.6,2X,F16.6,2X,F16.6)
    end do

    close(unit=9)

end subroutine writefile




end module
