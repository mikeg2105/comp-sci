!------------------------------------------------------------------------------
!
!  Program name:
!
!    inquire_iolength
!
!  Purpose:
!
!    This program shows the use of INQUIRE-by-output-list 
!    statement to determine the record length of an unformatted
!    output list. 
!
!------------------------------------------------------------------------------

PROGRAM inquire_iolength
  
  IMPLICIT NONE

  CHARACTER (LEN = 4) :: title = 'Mr.'
  CHARACTER (LEN = 15) :: name = 'Neil Smith'
  INTEGER :: age = 32
  CHARACTER (LEN = 35) :: address = '  2 Oak Lane, Manchester M23 4QD'
  INTEGER :: tel = 4455678  

  INTEGER :: rec_len		! record length to be inquired

  PRINT *, title, name, age, address, tel

  ! Assign the inquired record length to rec_len
  INQUIRE (IOLENGTH = rec_len) title, name, age, address, tel
  PRINT *, "The inquired record length = ", rec_len

  ! Open file for unformmated output using the inquired record length 
  OPEN(UNIT = 1, FILE = 'temp.dat', RECL = rec_len, FORM = 'UNFORMATTED')
  WRITE(1) title, name, age, address, tel
  CLOSE(UNIT = 1)  

  ! Open file for unformmated input
  OPEN(UNIT = 1, FILE = 'temp.dat', FORM = 'UNFORMATTED')
  READ(1) title, name, age, address, tel
  CLOSE(UNIT = 1)

  ! Print the list again (should be same as the first print)"
  PRINT *, title, name, age, address, tel

END PROGRAM
