program hello_user
    implicit none

    character(len=20) :: file
    character(len=20) :: name

    call get_file_name
    call read_name(trim(file), name)
    call write_greeting(trim(name))
contains
    subroutine get_file_name
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
    end subroutine

    subroutine read_name(file_name, name_read)
        character(len=*), intent(in) :: file_name
        character(len=*), intent(out) :: name_read

        integer :: unit

        open(newunit=unit, file=file_name, status="old")
        read(unit, '(A)') name_read
        close(unit)
    end subroutine

    subroutine write_greeting(name)
        character(len=*), intent(in) :: name

        integer :: unit

        open(newunit=unit, file="hello.txt", status="replace")
        write(unit, '(A)') trim(make_greeting(name))
        close(unit)
    end subroutine

    function make_greeting(name) result(greeting)
        character(len=*), intent(in) :: name
        character(len=100) :: greeting

        greeting = "Hello, " // name // "!"
    end function
end program
