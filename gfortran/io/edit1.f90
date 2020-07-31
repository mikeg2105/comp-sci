!------------------------------------------------------------------------------
!
!  Program name:
!
!    e_en_es_g_compare
!
!  Purpose:
!
!    To compare the differences between E, EN, ES and G
!    edit descriptors for output
!
!  Output result:
!
!    E:      0.123E+01    -0.500E+00     0.678E-02     0.988E+05
!    EN:     1.234E+00  -500.000E-03     6.780E-03    98.765E+03
!    ES:     1.234E+00    -5.000E-01     6.780E-03     9.877E+04
!    G:       1.23        -0.500         0.678E-02     0.988E+05
!
!------------------------------------------------------------------------------

PROGRAM e_en_es_g_compare

  IMPLICIT NONE
  REAL, DIMENSION(4) :: x = (/1.234, -0.5, 0.00678, 98765.4/)

  ! print array x
  PRINT '("E: ", 4E14.3)', x	! the 1st row is E  output
  PRINT '("EN:", 4EN14.3)', x	! the 2nd row is EN output
  PRINT '("ES:", 4ES14.3)', x	! the 3rd row is ES output
  PRINT '("G: ", 4G14.3)', x	! the 4th row is G  output

END PROGRAM e_en_es_g_compare
