program variables_example
    implicit none

    integer :: x
    integer :: y
    integer :: z

    print *, "x = ", x, ", y = ", y, ", z = ", z

    x = 3

    print *, "x = ", x, ", y = ", y, ", z = ", z

    y = 4

    print *, "x = ", x, ", y = ", y, ", z = ", z

    z = x + y

    print *, "x = ", x, ", y = ", y, ", z = ", z

    x = 1

    print *, "x = ", x, ", y = ", y, ", z = ", z
end program
