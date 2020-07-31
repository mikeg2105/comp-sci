!------------------------------------------------------------------------------
!
!  Program name:
!
!    munro_basis
!
!  Purpose
!
!    This program generates pseudo random real number between 0 and 4404, 
!    which provides the basis for you to write a complete program munro
!    to select all values greater than 3000 and find the number of them,
!    the maximum, the minimum and the average. 
!    
!------------------------------------------------------------------------------

PROGRAM munro_basis

  IMPLICIT NONE

  REAL, DIMENSION (5000) :: hills
  LOGICAL, DIMENSION (5000) :: bighills
  INTEGER, DIMENSION (1) :: iseed
  INTEGER :: nmunros

  REAL, ALLOCATABLE, DIMENSION (:) :: munros

  ! Generate pseudo random numbers.          

  iseed = 1
  CALL RANDOM_SEED (PUT=iseed)
  CALL RANDOM_NUMBER (hills(:))
  hills(:) = 4404.0*hills(:)

  ! Find number of Munros and allocate array munros:

  ! .......

  ! Pack Munros into array munros:

  ! .......

  ! Print number of, highest, lowest and average heights of Munros:

  ! ......

END PROGRAM munro_basis  
