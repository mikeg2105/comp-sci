program maths

    use prec, only:print_kind_info
    implicit none

    real :: r
    integer :: irad
    real, parameter :: pi=4.0*atan(1.0)

    r=3.0
    irad=4

    !area=pi*r**2
    print*, 'pi is',pi
    !print*, 'the area of a circle of radius ', r, 'is ',area(r)

    call print_kind_info()

end program maths
