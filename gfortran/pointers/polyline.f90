!------------------------------------------------------------------------------
!
!  Program name:
!
!    run_polyline
!
!  Purpose:
!
!    This program uses the module polyline to build and transverse a 
!    polyline linked list containing connected vertices.
!
!------------------------------------------------------------------------------

PROGRAM run_polyline

  USE polyline   	! makes all module definitions available

  IMPLICIT NONE

  REAL :: yval		! the y coordinate of a point
  TYPE(point) :: pt
  TYPE(vertex), POINTER :: current, first
  INTEGER :: status

  ! open files to read points and y value   

  OPEN (1, FILE = 'polyline.dat')

  ! build up a polyline list
   
  NULLIFY(first)		! initially nullify list (empty)

  DO
    READ(1, *) pt		! read point from channel 1 until 999.999 
    IF (ABS(pt%x - 999.999) < 0.001) EXIT
    ALLOCATE(current, STAT = status)		! create new node
    IF (status > 0) STOP 'Fail to allocate a new vertex'
    current%pt = pt     	! giving the point
    current%link => first	! point to previous one
    first => current		! update head of polyline list
  END DO
  
  WRITE(*, *) 'Printing the list for the first time'
  CALL print_list(first)
   
  WRITE(*, *)
  WRITE(*, *) 'Input the y value of the point to be deleted from the list'
  READ (*, *) yval		

  WRITE(*, *)
  WRITE(*, *) 'Looking up the list and deleting the node with the y value'
  CALL delentry(yval, first)

  WRITE(*, *)
  WRITE(*, *) 'Here is the new list:'
  CALL print_list(first)

END PROGRAM run_polyline
  


