program if_constructs
    implicit none

    if (.true.) print *, "This will be printed"
    if (.false.) print *, "This won't be printed"

    if (.true.) then
        print *, "This will also be printed"
    end if

    if (.false.) then
        print *, "This also won't be printed"
    end if

    if (.true.) then
        print *, "This will be printed too"
    else
        print *, "But this won't be"
    end if

    if (.false.) then
        print *, "This won't be printed"
    else
        print *, "But this will be"
    end if

    if (.false.) then
        print *, "This won't be printed"
    else if (.true.) then
        print *, "This one will be printed"
    else
        print *, "And this one won't be"
    end if
end program
