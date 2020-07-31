!---------------------------------------------------------------------
!
!  Program name:
!
!    i_b_o_z_compare
!
!  Purpose:
!
!    To compare the differences between I, B, O and Z
!    edit descriptors for integer output
!
!  Output result:
!
!  I:      1         8        15       245
!  B:      1      1000      1111  11110101
!  O:      1        10        17       365
!  Z:      1         8         F        F5
!
!  I:                                                              -15
!  B: 1111111111111111111111111111111111111111111111111111111111110001
!  O:                                           1777777777777777777761
!  Z:                                                 FFFFFFFFFFFFFFF1
!
!---------------------------------------------------------------------

PROGRAM i_b_o_z_compare

  IMPLICIT NONE
  INTEGER, DIMENSION(4) :: x = (/1, 8, 15, 245/)
  INTEGER :: i = -15

  ! print array x, 1st row is I output
  !                2nd row is B output
  !                3rd row is O output
  !                4th row is Z output

  PRINT '("I:", 4I10/"B:", 4B10/"O:", 4O10/"Z:", 4Z10)', x, x, x, x

  ! print i, note that it is negative number and 
  ! EL98 represent an integer with 64 bits by default.
  ! The correspondence of rows is same as the above.  

  PRINT '(/"I: ", I64/"B: ", B64/"O: ", O64/"Z: ", Z64)', i, i, i, i

END PROGRAM i_b_o_z_compare

