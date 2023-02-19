program while_loop
    implicit none

    integer :: i

    i = 100
    do while (i > 0)
        print *, "Can you guess how many times this will be printed?"
        i = i / 2
    end do
end program
