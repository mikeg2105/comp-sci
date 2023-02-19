program hello_user
    implicit none

    character(len=20) :: file
    character(len=20) :: name
    integer :: unit

    select case (command_argument_count())
    case (0)
        file = "name.txt"
    case (1)
        call get_command_argument(1, file)
    case default
        print *, "Too many command line arguments!"
        print *, "Usage: hello_user [file]"
        print *, "  where file is optional, and defaults to 'name.txt'"
        stop 1
    end select

    open(newunit=unit, file=trim(file), status="old")
    read(unit, '(A)') name
    close(unit)

    open(newunit=unit, file="hello.txt", status="replace")
    write(unit, '(A)') "Hello, " // trim(name) // "!"
    close(unit)
end program
