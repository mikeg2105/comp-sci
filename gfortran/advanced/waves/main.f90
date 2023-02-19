program main

    use io, only: write_file, read_file
    use prec, only:print_kind_info
    implicit none

    real :: r
    integer :: irad
    real, parameter :: pi=4.0*atan(1.0)
    character(len=20) :: filename
    real :: x1(4), y1(4)

    filename='results.out'
    x1 = (/1.,2.,3.,4./)
    y1 = x1**2
    r=3.0
    irad=4

    !area=pi*r**2
    print*, 'pi is',pi
    !print*, 'the area of a circle of radius ', r, 'is ',area(r)

    !call print_kind_info()
    !call write_file(filename,x1,y1)


    print*,'reinit data'
    x1=0.
    y1=0.

    print*,'x=',x1
    print*,'y=',y1

    print*,'now read file'

    call read_file(filename,x1,y1)

    print*,'x=',x1
    print*,'y=',y1



end program main
