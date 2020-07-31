!------------------------------------------------------------------------------
!
!  Program name:
!
!    assign_matrix
!
!  Purpose:
!
!    This program declares a 2-dimensional integer array, with extents (n,n), 
!    where n is set to 9 in a PARAMETER statement.
!
!    It uses DO loops to assign elements of the array to have values r c, 
!    where r is the row number and c is the column number, e.g., a(3,2) = 32, 
!    a(5,7) = 57. Then it writes the resulting matrix to the data file
!    matrix.dat for later use.  
!
!------------------------------------------------------------------------------

PROGRAM assign_matrix

  IMPLICIT NONE
  INTEGER, PARAMETER :: n = 9
  INTEGER, DIMENSION (n, n) :: ia
  INTEGER :: i, j

  ! Assign matrix ia with DO loops:

  WRITE(*, '("Assign with DO loops:"/)')

  DO i = 1, n
    DO j = 1, n
      ia(j,i) = i + j * 10
    END DO
  END DO

  WRITE(*,'(9I5)') (ia(i, :), i = 1, n)

  OPEN (UNIT = 1, FILE = 'matrix.dat')
  WRITE (1, *) ia
  CLOSE (UNIT = 1)
  
END PROGRAM assign_matrix
