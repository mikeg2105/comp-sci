program laxwend

use lparamdata
use laxsolver
use io, only: writeoutfile, read_file

implicit none
! use the lax-wendrof solution method to propaget a wave and a shock



include 'incparams.inc'




!This test is an MHD shocktube, where the right and left states are initialized to different values.
!The left state is initialized as (ρ, vx, vy, vz, By, Bz, p) = [1,0,0,1,0,1] and
!the right state [0.125,0,0,-1,0,0.1]. Bx=0.75 and γ = 2.
!The hydrodynamic portion of the initial conditions are the same as for the Sod shock tube problem.

filename='results.txt'

dx = xmax/(ni-1)



! 	Define the wavespeed
wavespeed = sqrt(fgamma*p0/rho0)



! 	Define time-domain
dt = courant*dx/wavespeed
tmax = 1.0d0   !15

nt=int(tmax/dt)
!nt=tmax/0.00001d0
!courant = wavespeed*dt/dx
print *, "courant, wavespeed, dt, nt"
print *, courant, wavespeed, dt, nt


call initpatarrays(ni,nt)
call initsolverarrays(ni,nt)




it=1




!initialise the configuration
do i=1,ni
    x(i)=(i-1)*dx
    if(i>76)then
        mx(i,1)=0
        rho(i,1)=rho0
        rhohm(i)=rho0
        rhohp(i)=rho0
        p(i,1)=p0
    else
        mx(i,1)=0
        rho(i,1)=rho0/8.0d0
        rhohm(i)=rho0/8.0d0
        rhohp(i)=rho0/8.0d0


        p(i,1)=p0/10.0d0
    endif
    e(i,1)=p(i,1)/((fgamma-1)*rho(i,1))
end do




do it=1,nt-1

    print *, "iteration: ",it
    !mx:1, rho:2, e:3
    ! for each field in turn
    !call flux()
    call flux(2, ni,it, mx(:,it), rho(:,it), e(:,it), p(:,it), frho, fgamma )
    call flux(1, ni,it, mx(:,it), rho(:,it), e(:,it), p(:,it), fmx, fgamma )
    call flux(3, ni,it, mx(:,it), rho(:,it), e(:,it), p(:,it), fe, fgamma )


    !for each field
    !call halfstep(timestep, deltax, denscurrent, dens, densh)
    !
    !    subroutine halfstep(timestep, deltax, denscurrent, dens, densh)
    !    subroutine halfstep(ni, nt, it, deltat, deltax, flux, field, fieldhdp,fieldhdm)

    call halfstep(ni, nt, it, dt, dx, fmx, mx, mxhm, mxhp)
    call halfstep(ni, nt, it, dt, dx, frho, rho, rhohm, rhohp)
    call halfstep(ni, nt, it, dt, dx, fe, e, ehm, ehp)

    ! for each field
    ! advance the final step
    !call advance()
    ! advance the density

    !    subroutine sadvance(ifield, ni, nt, it, deltat, deltax, field,  rhohm, rhohp, mxhm, mxhp, em, ep, fgamma)
    !mx:1, rho:2, e:3
    call sadvance(1, ni, nt, it, dt, dx, mx, rhohm, rhohp, mxhm, mxhp, ehm, ehp, fgamma)
    call sadvance(2, ni, nt, it, dt, dx, rho, rhohm, rhohp, mxhm, mxhp, ehm, ehp, fgamma)
    call sadvance(3, ni, nt, it, dt, dx, e, rhohm, rhohp, mxhm, mxhp, ehm, ehp, fgamma)

    mmx=sum(mx(:,it))/ni
    mrho=sum(rho(:,it))/ni
    me=sum(e(:,it))/ni
    p(:,it)=(fgamma-1)*rho(:,it)*e(:,it)
    mp=sum(p(:,it))/ni

    ! 	Define the wavespeed
    wavespeed = sqrt(fgamma*mp/mrho)
    ! 	Define time-domain
    dt = courant*dx/wavespeed
    dt=0.00001d0


    print*, dt, mmx,mrho,me
enddo ! end of time iterations

call writeoutfile(filename,x,mx,ni,nt)


contains

 !   function mean(it,ni,nt,field)
  !
  !
  !  endfunction






end program
