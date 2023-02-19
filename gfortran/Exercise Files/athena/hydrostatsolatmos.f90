program hydrostatsolatmos


    ! Routine uses tanh based solar temperature profile
    ! and hrdrostatic pressure balance

    !        height(m)        temp(K)     dens(kg)          pres()
    !        5955555.6       1599238.9   5.5055286e-12      0.12194448
    use statbalancemod, only: temp, writefile, hydropres, dens

    implicit none



    real, parameter :: deltah = 45454.545455
    real, parameter :: starth = 200.0
    integer, parameter :: npoints=132
    !real, parameter :: pi=4.0*atan(1.0)
    real,allocatable :: sh(:),sdens(:),spres(:),stemp(:)
    real :: hcurrent
    integer :: i


    allocate(sh(1:npoints))
    allocate(sdens(1:npoints))
    allocate(spres(1:npoints))
    allocate(stemp(1:npoints))

    ! initialise array values
    sh(:)=0
    sdens(:)=0
    spres(:)=0
    stemp(:)=0

    !perform the computation loop
    !hcurrent=deltah
    hcurrent=starth

    do i=1,npoints
        !calculate height
        sh(i)=hcurrent+deltah
        hcurrent=sh(i)

        ! calculate temperature
        stemp(i)=temp(sh(i))

    enddo


    hcurrent=starth
    do i=1,npoints
        !calculate height
        !sh(i)=hcurrent+deltah
        hcurrent=sh(i)


        !use hydrostatic balance to calculate  density
        spres(i)=hydropres(sh, i, npoints, deltah)
        print *, i,' height',sh(i),' dens ',spres(i)

        !compute pressure

        !write to ascii output file

!        print *, 'result'
    enddo


    do i=1,npoints

        sdens(i)=dens(sh(i),spres(i))

    enddo


    call writefile(sh,sdens,spres,stemp,npoints)

    print *, 'complete'
end program  hydrostatsolatmos
