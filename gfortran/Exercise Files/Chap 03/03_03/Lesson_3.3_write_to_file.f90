program hello_user
    implicit none

    character(len=20) :: name
    integer :: unit

    select case (command_argument_count())
    case (0)
        print *, "What's your name?"
        read(*,'(A)') name
    case (1)
        call get_command_argument(1, name)
    case default
        print *, "Too many command line arguments!"
        print *, "Usage: hello_user [name]"
        stop 1
    end select

    open(newunit=unit, file="hello.txt", status="replace")
    write(unit, '(A)') "Hello, " // trim(name) // "!"
    close(unit)
end program
