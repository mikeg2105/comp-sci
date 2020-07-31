!------------------------------------------------------------------------------
!
!  Program name:
!
!    vehicle
!
!  Purpose:
!
!    This program defines a derived type vreg, a vehicle registration 
!    number.
!
!------------------------------------------------------------------------------

PROGRAM vehicle

  IMPLICIT NONE
  TYPE vreg                      ! definition of vreg
    CHARACTER (LEN = 1) year      
    INTEGER             number
    CHARACTER (LEN = 3) place
  END TYPE vreg

  TYPE(vreg) mycar               ! mycar typed as vreg

  mycar = vreg('D', 608, 'PVR')  ! mycar assigned

  WRITE (*, *) 'printing by components'
  WRITE (*, '(1X, A1, I3, A3)') mycar%year, mycar%number, mycar%place
  WRITE (*, *) 'printing by structure'
  WRITE (*, *) mycar

END PROGRAM vehicle
