!------------------------------------------------------------------------------
!
!  Program name:
!
!    model_numbers
!
!  Purpose:
!
!    This program shows the use of intrinsics for floating point 
!    manipulation and inquiry. 
!
!------------------------------------------------------------------------------

PROGRAM model_numbers

  IMPLICIT NONE

  INTEGER, PARAMETER :: ik = SELECTED_REAL_KIND(10)
  REAL(ik), PARAMETER, DIMENSION(2) :: a = (/ 1.0_ik, 3.14159265358979323846_ik /)
  REAL(ik) :: pi, one = a(1)
  INTEGER :: i, n = 3

  DO i = 1, 2
    pi = a(i)

    ! Using Intrinsics for floating point manipulation:

    WRITE (*,*)
    WRITE (*,*) 'Original value: ', pi
    WRITE (*,*)
    WRITE (*,*) 'Exponent      : ', EXPONENT (pi)
    WRITE (*,*) 'Fraction      : ', FRACTION (pi)
    WRITE (*,*) 'Nearest above : ', NEAREST (pi,one)
    WRITE (*,*) 'Rrspacing     : ', RRSPACING (pi)
    WRITE (*,*) 'Scale         : ', SCALE (pi,n)
    WRITE (*,*) 'Set_exponent  : ', SET_EXPONENT (pi,n)
    WRITE (*,*) 'Spacing       : ', SPACING (pi)
    WRITE (*,*)

    ! Using numeric enquiry intrinsics.

    WRITE (*,*) 'Digits        : ', DIGITS(pi)
    WRITE (*,*) 'Epsilon       : ', EPSILON(pi)
    WRITE (*,*) 'Huge          : ', HUGE(pi)
    WRITE (*,*) 'Max exponent  : ', MAXEXPONENT(pi)
    WRITE (*,*) 'Min exponent  : ', MINEXPONENT(pi)
    WRITE (*,*) 'Precision     : ', PRECISION(pi)
    WRITE (*,*) 'Radix         : ', RADIX(pi)
    WRITE (*,*) 'Tiny          : ', TINY(pi)
    WRITE (*,*)

END DO

END PROGRAM model_numbers





























