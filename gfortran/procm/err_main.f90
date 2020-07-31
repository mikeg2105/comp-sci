!------------------------------------------------------------------------------
!
!  Program name:
!
!    errors
!
!  Purpose:
!
!    Look at this program.  What is wrong with it?  Rewrite it in a
!    better way in Fortran 90.
!
!------------------------------------------------------------------------------


PROGRAM errors

! errors.f90          - finds ratio of two values

  IMPLICIT NONE

  INTEGER, PARAMETER :: iki = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER :: ikr = SELECTED_REAL_KIND(10,200)

  REAL(ikr) :: ia, ib, ic, id
  INTEGER(iki) :: a, b, c, d

  a = 3_ikr
  b = 2_ikr    

  ia = 3_iki
  ib = 2_iki

! Find ratios:

  CALL ratio(a, b, c, ia, ib, ic)
  CALL ratio(b, a, d, ib, ia, id)

! Multiply ratios to get unity:

  WRITE(*,*)  c*d, ic*id

END PROGRAM errors
  
