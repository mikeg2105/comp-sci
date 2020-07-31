!------------------------------------------------------------------------------
!
!  Program name:
!
!    kind_int
!
!  Purpose:
!
!    This program uses SELECTED_INT_KIND to find the kind values for
!    integer variables on the system.
!
!------------------------------------------------------------------------------

PROGRAM kind_int

  IMPLICIT NONE
  INTEGER :: i

! SELECTED_INT_KIND(R)    returns kind value for integer in range
!                         -10**R < integer < 10**R
!
! SELECTED_REAL_KIND(P,R) returns kind value for real with P decimal
!                         places of precision and exponent range R.

  WRITE(*,*)
  WRITE(*,*) '     I      Kind Value' 
  WRITE(*,*) '            for int.  '
  WRITE(*,*) '            with range'
  WRITE(*,*) '            +/- 10**I '
  WRITE(*,*)

  DO i = 1, 10
    WRITE(*, '(I6, I12)') i, &
             SELECTED_INT_KIND(i)
  END DO

END PROGRAM kind_int
