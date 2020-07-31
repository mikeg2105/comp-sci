  PROGRAM stats
!
! program to calculate the mean and standard deviation of a set of data values.
! 
  IMPLICIT NONE
  INTEGER, PARAMETER :: NMAX=200 
  INTEGER :: n , i          
  REAL :: x(NMAX) , xmean , xstd , sum
  
  print * , 'Enter the number of data values to be read'  
  read *, n
  print *, 'Enter the data values:' 
  read (*,*) ( x(i) , i = 1 , n )  
                                   
  
!
! calculate the mean
!  
   sum = 0.0              
   DO i = 1  , n
       sum = sum + x(i)   
   ENDDO
   xmean = sum / n  
!
! calculate the standard deviation
!                
   sum = 0.0
   DO i = 1 , n 
       sum = sum + (x(i)-xmean)**2
   ENDDO
   xstd = sqrt ( sum /(n-1 ) )  
   write (*,100) xmean , xstd 
! above is a formatted write statement.
! Statement labeled 100 is used for output.
! Test different effects by uncommenting one and only 
! one of the various versions of the Format statement 100 below.  
 100 format ( 'MEAN VALUE=' , F12.4 ,5X,'STANDARD DEVIATION=',F12.6 )
! 100 format ( 'Mean Value:',E12.4 / 'Standard Deviation:' ,E12.4 )
! 100 format ( 'Mean Value:',T20,G12.3 /'Standard Deviation:',T20, F5.0 )
 
  END
