!------------------------------------------------------------------------------
!
!  Program name:
!
!    convert
!
!  Purpose:
!
!    This program shows the use of elemental numeric intrinsics for 
!    converting real in different ways (AINT, ANINT, INT, NINT, 
!    CELLING, FLOOR).
!
!------------------------------------------------------------------------------

PROGRAM convert

  IMPLICIT NONE

  INTEGER, PARAMETER :: ik = SELECTED_REAL_KIND(10)
  REAL(ik) :: pi = -3.141592653589793_ik
  INTEGER :: ia

  WRITE(*,*) 'Kind value for reals with > 10 d.p.: ', ik
  WRITE(*,*)
  WRITE(*,*) 'Original value: ', pi
  WRITE(*,*)

  ! FORTRAN 77:

  WRITE(*,*) 'AINT    - Truncate to real    : ', AINT(pi)
  WRITE(*,*) 'ANINT   - Nearest real        : ', ANINT(pi)
  WRITE(*,*) 'INT     - Truncate to integer : ', INT(pi)
  WRITE(*,*) 'NINT    - Nearest integer     : ', NINT(pi)

  ! Fortran 90:

  WRITE(*,*) 'CEILING - Next larger integer : ', CEILING(pi)
  WRITE(*,*) 'FLOOR   - Next smaller integer: ', FLOOR(pi)

END PROGRAM convert
