!module containing geometry functions

module geometry

    implicit none
    real, parameter :: pi=4.*atan(1.0)
    public :: area, pi

    private


    contains

    !
    ! a function to calculate the area of a circle
    !
    real function area(r)
    real, intent(in) :: r

    area = pi*r**2
    end function area

end module geometry
