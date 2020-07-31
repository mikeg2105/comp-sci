!------------------------------------------------------------------------------
!
!  Program name:
!
!    non_advancing_io
!
!  Purpose:
!
!    This program shows the use of non advancing I/O to suppress 
!    line feed after the input prompt.
!
!    This program also shows the use of an allocatable array and 
!    array intrinsic function SPREAD.
!
!------------------------------------------------------------------------------

PROGRAM non_advancing_io

  IMPLICIT NONE
  INTEGER :: n, i
  INTEGER, DIMENSION (:,:), ALLOCATABLE :: ia

  ! suppress line feed after input prompt by non-advancing io
  WRITE(*, '("Extent of matrix? ")', ADVANCE = 'NO')

  READ(*, '(I4)') n		! read in rank of matrix

  ALLOCATE (ia(n, n))		! allocate array ia 

  ! SPREAD intrinsic function is used to construct an array
  ! by duplicating an array section
  ia = SPREAD ( (/ (i, i = 1, n) /), DIM = 1, NCOPIES = n)	&     
     + SPREAD ( (/ (i, i = 1, n) /), DIM = 2, NCOPIES = n)

  ! output the values of the array ia
  DO i = 1, n
    WRITE(*, '(10I5)') ia(:, i)
  END DO

END PROGRAM non_advancing_io
