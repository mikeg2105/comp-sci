!------------------------------------------------------------------------------
!
!  Program name:
!
!    spread
!
!  Purpose:
!
!    This program shows how to use the array intrinsic function SPREAD to 
!    assign an array.
!
!------------------------------------------------------------------------------

PROGRAM spread1

  IMPLICIT NONE
  INTEGER, PARAMETER :: n = 5
  INTEGER :: i
  REAL, DIMENSION (n,n) :: ra

  WRITE(*,'("Using SPREAD with DIM = 2:"//)')
  ra = REAL( SPREAD( (/ (i, i = 1, n) /), DIM = 1, NCOPIES = n) ) 
  WRITE(*,'(5F8.2)') (ra(:,i), i = 1, n)

  WRITE(*,'("Using SPREAD with DIM = 1:"//)')
  ra = REAL( SPREAD( (/ (i, i = 1, n) /), DIM = 2, NCOPIES = n) ) 
  WRITE(*,'(5F8.2)') (ra(:,i), i = 1, n)

END PROGRAM spread1



