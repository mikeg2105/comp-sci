!------------------------------------------------------------------------------
!
!  Program name:
!
!    ch_elemental
!
!  Purpose:
!
!    This program shows the use of elemental intrinsics for 
!    characters (ADJUSTL, ADJUSTR, INDEX, LEN_TRIM, SCAN, VERIFY).
!
!------------------------------------------------------------------------------

PROGRAM ch_elemental

  IMPLICIT NONE
  CHARACTER (LEN = 2) :: chstars = '**'
  CHARACTER (LEN = 40) :: &
    char = '     5 Blanks,  30 chars,  5 blanks     '

  WRITE (*,*) 'Original : ', chstars, char, chstars
  WRITE (*,*) 'Adjustl  : ', chstars, ADJUSTL (char), chstars
  WRITE (*,*) 'Adjustr  : ', chstars, ADJUSTR (char), chstars
  WRITE (*,*) 'Index    : ', INDEX(char,'30')
  WRITE (*,*) 'Len_trim : ', LEN_TRIM (char)
  WRITE (*,*) 'Scan     : ', SCAN (char, 'cha')
  WRITE (*,*) 'Verify   : ', VERIFY (char, 'blank')

END PROGRAM ch_elemental
