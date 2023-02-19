program squares
    implicit none

    integer :: i
    integer :: the_squares(10)

    the_squares = 42
    print *, the_squares

    do i = 1, 10
        the_squares(i) = i**2
    end do

    do i = 10, 1, -1
        print *, the_squares(i)
    end do

    print *, the_squares(2:4)
end program
