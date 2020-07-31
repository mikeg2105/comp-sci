!------------------------------------------------------------------------------
!
!  Program name:
!
!    simple_linked_list
!
!  Purpose:
!
!    This program demonstrates how to build and transverse a 
!    simple linked list by using pointers and dynamic allocation
!    facilities.
!
!    Notice the linked list stores the typed-in numbers in 
!    reverse order.
!
!------------------------------------------------------------------------------

PROGRAM simple_linked_list

  IMPLICIT NONE
  TYPE node
    INTEGER :: value			! data filed
    TYPE (node), POINTER :: next	! pointer field
  END TYPE node
  INTEGER :: num, status
  TYPE (node), POINTER :: list, current

  ! build up the list

  NULLIFY(list)				! initially nullify list (empty)

  PRINT *, 'Type-in an integer number to build a linked list (0 to terminate)'

  DO
    READ *, num				! read num from keyboard 
    IF (num == 0) EXIT			! until 0 is entered
    ALLOCATE(current, STAT = status)    ! create new node
    IF (status > 0) STOP 'Fail to allocate a new node'
    current%value = num			! give the value
    current%next => list		! point to previous one
    list => current			! update head of list
  END DO

  ! transverse the list and print the values

  PRINT *, 'Transverse the list built up and print the values'

  current => list			! make current as alias of list
  DO
    IF (.NOT. ASSOCIATED(current)) EXIT	! exit if null pointer
    PRINT *, current%value		! print the value
    current => current%next 		! make current alias of next node
  END DO

END PROGRAM simple_linked_list




