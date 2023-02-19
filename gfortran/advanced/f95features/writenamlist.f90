! example.f90
program main
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
    implicit none

    type :: person_type
        integer           :: id
        character(len=32) :: name
        integer           :: age
    end type person_type

    integer           :: x, y
    real              :: r(2)
    type(person_type) :: person

    ! Set initial values.
    person = person_type(-1, 'Jane Does', 0)
    x = 0; y = 1
    r = [ 2.0, 0.0 ]

    ! Read from file.
    call write_namelist('namelist.nml', person, x, y, r)

    ! Output some values.
    print '("PERSON ID:   ", i0)', person%id
    print '("PERSON NAME: ", a)',  person%name
    print '("PERSON AGE:  ", i0)', person%age
contains
    subroutine write_namelist(file_path, person, x, y, r)
        !! Reads Namelist from given file.
        character(len=*),  intent(in)    :: file_path
        type(person_type), intent(inout) :: person
        integer,           intent(inout) :: x, y
        real,              intent(inout) :: r(2)
        integer                          :: fu, rc

        ! Namelist definition.
        namelist /EXAMPLE/ x, y, r, person

        ! Check whether file exists.
        inquire (file=file_path, iostat=rc)

        if (rc /= 0) then
            write (stderr, '("Error: input file ", a, " does not exist")') file_path
            return
        end if

        ! Open and read Namelist file.
        open (action='write', file=file_path, newunit=fu)
        write (nml=EXAMPLE, unit=fu)
        
        !if (rc /= 0) write (stderr, '("Error: invalid Namelist format")')

        close (fu)
    end subroutine write_namelist
end program main
