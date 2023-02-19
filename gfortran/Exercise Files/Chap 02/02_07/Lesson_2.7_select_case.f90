program select_case
    implicit none

    select case (1+2)
    case (2+1, 2+2)
        print *, "was 3 or 4"
    end select

    select case ("Hi")
    case ("Hello")
        print *, "Hello"
    case default
        print *, "Default"
    end select

    select case (3)
    case (:4)
        print *, "Less than or equal to 4"
    end select
end program
