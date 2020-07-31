!------------------------------------------------------------------------------
!
!  Program name:
!
!    circledef4
!
!  Purpose:
!   
!    This program defines the derived types of point and circle,  
!    and shows how to assign a value to a derived type directly.
!
!------------------------------------------------------------------------------

PROGRAM circledef4 

  IMPLICIT NONE
  TYPE point         ! defines type point
    REAL x,y
  END TYPE point

  TYPE circle        ! defines type circle
    TYPE(point) centre
    REAL        radius
  END TYPE circle

  TYPE(circle) one   ! one is type circle

  ! Next statement assigns value to one (of derived type circle)

  one = circle(point(10.,20.),5.) 

  PRINT '(//)'
  PRINT *,'Print circle by derived type'
  WRITE (*,*) one
  PRINT '(//)'

END PROGRAM circledef4

