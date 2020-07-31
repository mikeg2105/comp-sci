!------------------------------------------------------------------------------
!
!  Program name:
!
!    namelist
!
!  Purpose:
!
!    This program shows the use of NAMELIST for input and output
!
!  Note:
!
!    This program reads the data file namelist.dat.
!
!    In OPEN statements, the DELIM specifier is used, which is
!    necessary to delimit any character strings so that the
!    output sent to a file could be read subsequently by another 
!    program, or even by the same one, using directed or namelist
!    formatting.         
!
!------------------------------------------------------------------------------

PROGRAM nmlist

  IMPLICIT NONE
  INTEGER :: quantity
  CHARACTER(LEN=6) :: colour(4)=(/ 'red   ', 'black ', 'green ', 'blue  ' /)
  CHARACTER(LEN=6) :: size(3)=(/ ' small', 'medium', ' large' /)
  NAMELIST /clothes/ quantity, colour, size	! define namelist group

  OPEN(6, DELIM = "APOSTROPHE")
  WRITE(6, NML=clothes)	! use namelist name as format for output

  OPEN(1, FILE='namelist.dat', DELIM = "APOSTROPHE")
  READ(1, NML=clothes)	! use namelist name as format for input
  CLOSE(1)

  WRITE(6, NML=clothes)	! output the newly read namelist 
  CLOSE(6)

END PROGRAM nmlist
