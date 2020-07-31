!------------------------------------------------------------------------------
!
!  Program name:
!
!    noconvert
!
!  Purpose:
!
!    This program shows the use of elemental numeric intrinsics that 
!    do not convert (MOD, MODULO)
!
!------------------------------------------------------------------------------


PROGRAM noconvert

  IMPLICIT NONE

  INTEGER :: a = -100, p = 7

  WRITE(*,*) 'Original values of a and p: ', a, p
  WRITE(*,*)

  ! FORTRAN 77:

  WRITE(*,*) 'MOD     - Remainder of a modulo p: ', MOD(a,p)

  ! Fortran 90:

  WRITE(*,*) 'MODULO  - a modulo p             : ', MODULO(a,p)

END PROGRAM noconvert
