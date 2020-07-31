!------------------------------------------------------------------------------
!
!  Program name:
!
!    convert_case
!
!  Purpose:
!
!    This program shows how to use DO and CASE control constructs. 
!    It reads in a character string, changes the case of letters and 
!    writes out new string.
!
!  Note:
!
!    Need to know the difference in the collation sequence of the upper 
!    and lower case characters - use position of upper and lower case A: 
!
!                        IACHAR('A') - IACHAR('a')
!
!------------------------------------------------------------------------------

PROGRAM convert_case

  IMPLICIT NONE

  CHARACTER (LEN = 80) :: string
  INTEGER :: i, upper_to_lower, len_string

  upper_to_lower = IACHAR("a") - IACHAR("A")

  WRITE (*,*) 'Please enter string of up to 80 characters'

  ! Read in string (up to 80 characters):
  READ (*,'(A)') string

  ! Find length of string excluding trailing blanks
  len_string = LEN_TRIM(string)

  DO i=1,len_string

  ! Use CASE construct to change case of the ith character.
  ! Use two cases, one for upper to lower case and one for 
  ! lower to upper case. 
 
    !!
 
    REMOVE THIS LINE AND FILL IN CODE HERE

    !!

  END DO

  ! Write out converted string
  WRITE (*,'(A)') string

END program convert_case















