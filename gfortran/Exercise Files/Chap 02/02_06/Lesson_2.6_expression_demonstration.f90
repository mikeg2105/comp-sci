program expression_demonstration
    implicit none

    integer, parameter :: x = 3
    integer :: y
    integer :: z

    y = 5
    z = (sin(0.5 - 0.5)*x + fib(5)*y - 2**4)/5
    print *, z
contains
    pure recursive function fib(n) result(fib_)
        integer, intent(in) :: n
        integer :: fib_

        if (n == 1 .or. n == 2) then
            fib_ = 1
        else
            fib_ = fib(n-1) + fib(n-2)
        end if
    end function
end program
