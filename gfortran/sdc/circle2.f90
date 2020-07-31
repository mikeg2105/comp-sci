!------------------------------------------------------------------------------
!
!  Program name:
!
!    circledef2
!
!  Purpose:
!
!    This program derives a defined type for a circle and
!    reads in the values by components from a data file
!    opened on channel 1.
!
!------------------------------------------------------------------------------

PROGRAM circledef2  
                    
  IMPLICIT NONE      ! no implicit typing

  TYPE point         ! defines type point
    REAL x,y
  END TYPE point

  TYPE circle        ! defines type circle
    TYPE(point) centre
    REAL        radius
  END TYPE circle

  TYPE(point) o      ! o is type point
  TYPE(circle) one   ! one is type circle

  OPEN (1,FILE='circle.dat')  ! opens file on channel 1

  READ (1,'(2F6.2)') o%x,o%y     ! reads in point o by component
  one%centre = o

  READ (1,'(F6.2)') one%radius   ! reads in radius by component

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

END PROGRAM circledef2
