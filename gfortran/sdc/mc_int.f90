!------------------------------------------------------------------------------
!
!  Program name:
!
!    machine_constants_int
!
!  Purpose:
!
!    This program uses numeric intrinsic functions to find some machine 
!    constants for integer variables with different kind values on the system.
!
!------------------------------------------------------------------------------

PROGRAM machine_constants_int

  IMPLICIT NONE

  INTEGER, PARAMETER, DIMENSION(3) :: ikvl = (/ 1,2,3 /)      ! Cray CS6400 
  INTEGER (KIND = ikvl(1)) i1
  INTEGER (KIND = ikvl(2)) i2
  INTEGER (KIND = ikvl(3)) i3

  WRITE (*,*) &
    'Machine constants for integer variables with different kind values'
  WRITE (*,*)
  WRITE (*,*) 'Default interger kind value: ', KIND(1)

  WRITE (*,'(/''   ikv   Digits           Huge    Range''/)')
  WRITE(*,'(i5,i7,i20,i5)') ikvl(1), DIGITS(i1), HUGE(i1), RANGE(i1)
  WRITE(*,'(i5,i7,i20,i5)') ikvl(2), DIGITS(i2), HUGE(i2), RANGE(i2)
  WRITE(*,'(i5,i7,i20,i5)') ikvl(3), DIGITS(i3), HUGE(i3), RANGE(i3)

END PROGRAM machine_constants_int
