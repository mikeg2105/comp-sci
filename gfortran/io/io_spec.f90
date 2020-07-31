!------------------------------------------------------------------------------
!
!  Program name:
!
!    io_specifier
!
!  Purpose:
!
!    To demonstrate some of the new I/O facilities of Fortran 90:
!      o New statement specifiers: POSITION, ACTION, ADVANCE, SIZE
!      o Non-advancing I/O
!
!  Note:
!
!    This program will write to or read from the data file 
!    io_specifier.dat
!
!------------------------------------------------------------------------------

PROGRAM io_specifier

  IMPLICIT NONE
  INTEGER, PARAMETER :: n = 10
  INTEGER :: i, isize
  INTEGER, DIMENSION (n) :: ix

  OPEN (7, FILE = 'io_specifier.dat', & 
	   POSITION = 'ASIS',         &
	   ACTION = 'WRITE')

  ! Write 10 values to unit 7:
  WRITE(7, '(10I5)') (i, i = 1, n)

  CLOSE(7)

  ! Rewind the file
  OPEN (7, FILE = 'io_specifier.dat', &
	   POSITION = 'REWIND',       &
	   ACTION = 'READWRITE')

  ! Read first 5 elements of first record:
  READ(7, '(5I5)', ADVANCE = 'NO', SIZE = isize) (ix(i), i = 1, 5)
  PRINT *, "isize = ", isize, "  ix(1:10) = ", ix

  ! Read next 3 elements of first record:
  READ(7, '(3I5)', ADVANCE = 'NO', SIZE = isize) (ix(i), i = 6, 8)
  PRINT *, "isize = ", isize, "  ix(1:10) = ", ix

  ! Write another 10 elements to unit 7. File was positioned after
  ! the 8th value, so there will now be 18 values.
  WRITE(7, '(10I5)') (-i, i = 1, n)

  ! Implicit ADVANCE='YES' in WRITE, so file now positioned at 
  ! beginning of 2nd record. Rewind again:
  REWIND(7)

  ! Read first 10 elements of first record 
  ! (should be 1, 2, 3, 4, 5, 6, 7, 8, -1, -2)
  READ(7, '(10I5)', ADVANCE = 'NO', SIZE = isize) (ix(i), i = 1, 10)
  PRINT *, "isize = ", isize, "  ix(1:10) = ", ix

  CLOSE(7)
  
END PROGRAM io_specifier
