!------------------------------------------------------------------------------
!
!  Program name:
!
!    generic
!
!  Purpose:
!
!    This program demonstrates how to define a generic procedure for 
!    swapping reals, integers and logicals using external procedures.
!
!------------------------------------------------------------------------------

PROGRAM generic

  IMPLICIT NONE
  REAL :: a, b
  INTEGER :: ia, ib
  LOGICAL :: la, lb

  INTERFACE swap		! generic name
    SUBROUTINE swapreal(u,v)
      REAL, INTENT(INOUT) :: u, v
    END SUBROUTINE swapreal

    SUBROUTINE swapint(u,v)
      INTEGER, INTENT(INOUT) :: u, v
    END SUBROUTINE swapint

    SUBROUTINE swaplog(u,v)
      LOGICAL, INTENT(INOUT) :: u, v
    END SUBROUTINE swaplog
  END INTERFACE

  ! swap two reals
  a = 1.0
  b = 2.0
  WRITE(*,*) a, b
  CALL swap(a,b)
  WRITE(*,*) a, b
  WRITE(*,*)

  ! swap two integers
  ia = 1
  ib = 2
  WRITE(*,*) ia, ib
  CALL swap(ia,ib)
  WRITE(*,*) ia, ib
  WRITE(*,*)

  ! swap two logicals
  la = .TRUE.
  lb = .FALSE.
  WRITE(*,*) la, lb
  CALL swap(la,lb)
  WRITE(*,*) la, lb
  WRITE(*,*)

END PROGRAM generic

SUBROUTINE swapreal(u,v)	! external procedure to swap two reals

  IMPLICIT NONE
  REAL, INTENT(INOUT) :: u, v
  REAL :: temp

  temp = u
  u = v
  v = temp

END SUBROUTINE swapreal

SUBROUTINE swapint(u,v)		! external procedure to swap two integers

  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: u, v
  INTEGER :: temp

  temp = u
  u = v
  v = temp

END SUBROUTINE swapint

SUBROUTINE swaplog(u,v)		! external procedure to swap two logicals

  IMPLICIT NONE
  LOGICAL, INTENT(INOUT) :: u, v
  LOGICAL :: temp

  temp = u
  u = v
  v = temp

END SUBROUTINE swaplog






