      PROGRAM STARTSHOOTING
!  THIS PROVIDES A SOLUTION FOR A BOUNDARY-VALUE PROBLEM FOR A
!  FIRST-ORDER ODE WITH ONE UNKNOWN PARAMETER TO BE DETERMINED.
!  THE ODE IS OF THE FORM     Dy/Dx = f(x,Q)  WHERE Q IS THE 
!  UNKNOWN PARAMETER.   THE BOUNDARY CONDITIONS ARE y = Ya FOR
!  x = a AND y = Yb FOR X = b.   
!  THE FUNCTION f(x,Q) IS PROVIDED AS A FUNCTION STATEMENT.
      REAL Y(0:1000),W(4),DY(4)
      CHARACTER ANS*1
      DATA W/0.0,0.5,0.5,1.0/
      FUN(X,Q)=-15.915494*Q/(2-X)**2
      WRITE(6,'('' INPUT THE FIRST BOUNDARY CONDITION AS a,Ya'')')
      WRITE(6,'('' WHERE y = Ya WHEN x = a.'')')
      READ(5,*)XI,Y(0)
      WRITE(6,'('' INPUT THE FINAL BOUNDARY CONDITION AS b,Yb'')')
      WRITE(6,'('' WHERE y = Yb WHEN x = b.'')')
      READ(5,*)XF,YF
      WRITE(6,'('' INPUT THE INTEGRATION STEP LENGTH h IN THE '')')
      WRITE(6,'('' FORM OF AN INTEGER N WHERE h = [b - a]/N. '')')
      WRITE(6,'('' IF YOU WANT OUTPUT OF Y AT INTERVALS OF [B-A]/M'')')
      WRITE(6,'('' THEN MAKE N A MULTIPLE OF M.'')')
      READ(5,*)N
      H=(XF-XI)/N
      WRITE(6,'('' INPUT ESTIMATE OF THE UNKNOWN PARAMETER, Q.'')')
      READ(5,*)Q
!  THE RUNGE-KUTTA INTEGRATION NOW BEGINS.F
      DO 1 I=1,N
      X=(I-1)*H
      DO 2 J=1,4
      XX=X+W(J)*H
      DY(J)=H*FUN(XX,Q)
    2 CONTINUE
      Y(I)=Y(I-1)+(DY(1)+DY(4)+2*(DY(2)+DY(3)))/6.0
    1 CONTINUE
      WRITE(6,100)Y(N)
  100 FORMAT(20H THE VALUE OF Yb IS ,F8.2)
    5 WRITE(6,'('' OUTPUT INTERMEDIATE VALUES OF Y? [Y/N]'')')
      READ(5,50)ANS
   50 FORMAT(A1)           
      IF(ANS.EQ.'Y'.OR.ANS.EQ.'y')GOTO 3
      IF(ANS.EQ.'N'.OR.ANS.EQ.'n')GOTO 4
      GOTO 5
    3 WRITE(6,'('' DO YOU WANT PRINTED OUTPUT? [Y/N]'')')
      READ(5,50)ANS
      IF(ANS.EQ.'Y'.OR.ANS.EQ.'y')THEN
      IOUT=9
      OPEN(UNIT=9, FILE='mydoc.txt')
      GOTO 7
      ENDIF
      IF(ANS.EQ.'N'.OR.ANS.EQ.'n')THEN
      IOUT=6
      GOTO 7
      ENDIF
      GOTO 3
    7 WRITE(6,'('' INPUT M WHERE OUTPUT INTERVALS ARE H*N/M '')')
      READ(5,*)M
      K=N/M
      DO 6 I=0,N,K
      XX=I*H
      WRITE(IOUT,200)XX,Y(I)
  200 FORMAT(F6.3,2X,F8.2)
    6 CONTINUE
    4 STOP
 END PROGRAM STARTSHOOTING

 FUNCTION myfunc1(x,q)	! external procedure to swap two reals

  IMPLICIT NONE
  REAL, INTENT (IN):: x,q
  REAL :: myfunc1
   myfunc1=  -15.915494*q/(2-x)**2

 
END FUNCTION myfunc1
