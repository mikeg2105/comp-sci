program hello_user
    implicit none

    character(len=20) :: name

    print *, "What's your name?"
    read(*,'(A)') name
    print *, "Hello, " // trim(name) // "!"
end program
