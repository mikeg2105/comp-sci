!------------------------------------------------------------------------------
!
!  Program name:
!
!    circledef1
!
!  Purpose:
!   
!    This program defines a derived type for a circle and
!    shows how to pass a user-defined type to an internal
!    subprogram.
!
!------------------------------------------------------------------------------

PROGRAM circledef1
                  
  IMPLICIT NONE
  TYPE point         ! defines type point
    REAL x,y
  END TYPE point

  TYPE circle        ! defines type circle
    TYPE(point) centre
    REAL        radius
  END TYPE circle

  TYPE(point) o             ! o is type point
  TYPE(circle) one,double   ! one and double are type circle

  o%x = 1.0          ! initialises by components to avoid confusion
  o%y = 2.0          ! with complex constant
  one%centre = o     ! initialises centre component of circle one
  one%radius = 0.5   ! initialises the radius of circle one

  PRINT '(/)'
  PRINT *,'Print centre of circle by components'
  WRITE (*,'(2F6.2,/)') o%x,o%y

  PRINT *,'Print centre of circle by derived type'
  WRITE (*,*) o

  PRINT '(//)'
  WRITE (*,*) 'Print circle by components'
  WRITE (*,'(3F6.2,/)') one

  PRINT *,'Print circle by derived type'
  WRITE (*,*) one
  PRINT '(//)'

  CALL circlex2(one,double)

  PRINT *,'Now follows the circle after being doubled in radius'
  PRINT '(//)'

  WRITE (*,*) 'Print circle by components'
  WRITE (*,'(3F6.2,/)') double

  PRINT *,'Print circle by derived type'
  WRITE (*,*) double
  PRINT '(/)'

CONTAINS

  SUBROUTINE circlex2(one,two)
    TYPE(CIRCLE) one,two
    two%centre = one%centre
    two%radius = one%radius * 2.0
  END SUBROUTINE circlex2 

END PROGRAM circledef1

