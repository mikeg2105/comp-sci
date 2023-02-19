program indexed_loops
    implicit none

    integer :: i

    do i = 1, 5
        print *, i
    end do

    do i = 2, 1
        print *, "This won't be printed even once"
    end do

    do i = 10, 1, -2
        print *, i
    end do

    do i = 1, 2, -1
        print *, "This won't be printed either"
    end do
end program
