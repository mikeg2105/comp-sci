  PROGRAM QES
!
! program to calculate the roots of a quadratic eqn.
! 
  IMPLICIT NONE 
  INTEGER :: a, b, c, D 
  REAL :: RealP, ImagP 
  a = 1
repeat:  do while ( a /= 0 )
    PRINT *, ' Type in three whole numbers for a, b and c or all three zeros to finish'

    READ *, a, b, c 
    IF (a /= 0) THEN 
! Calculate discriminant 
    D = b*b - 4*a*c
    IF (D == 0) THEN        ! one root 
       PRINT *,  'Root is ' ,  -b/(2.0*a) 
    ELSE IF (D > 0) THEN   ! real roots 
       PRINT *, 'Roots are' , (-b+SQRT(REAL(D)))/(2.0*a) , & 
              ' and ' ,  (-b-SQRT(  REAL(D) ) ) / (2.0*a) 
    ELSE ! complex roots 
       RealP = - b / (2.0*a) 
       ImagP = (SQRT(REAL(-D))/(2.0*a)) 
       PRINT *,   '1st Root' ,   RealP,    '+'   ,   ImagP,  'i' 
       PRINT *,   '2nd Root' ,   RealP,    '¡'    ,  ImagP,  'i' 
    ENDIF 
  ELSE                            ! a is equal to 0 
    PRINT*, 'Finished ...'
  ENDIF 
  end do repeat
  END PROGRAM QES
