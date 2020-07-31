  PROGRAM stats
!
! program to calculate the mean and standard deviation of a set of data values.
! 
  IMPLICIT NONE
  INTEGER, PARAMETER :: NMAX=200 
  INTEGER :: n , i          
  REAL :: x(NMAX) , xmean , xstd , sum
  
  
  call readdata
  call calcmean
  call calcstd
  write (*,100) xmean , xstd 
100 format ( 'MEAN VALUE=' , F12.4 ,5X,'STANDARD DEVIATION=',F12.6 )
  stop
  
  contains
    subroutine readdata
  
      print * , 'Enter the number of data values to be read'  
      read *, n
      print *, 'Enter the data values:' 
      read (*,*) ( x(i) , i = 1 , n )  
    end subroutine readdata
    
    subroutine calcmean                               
  
!
! calculate the mean
!     
       sum = 0.0              
       DO i = 1  , n
         sum = sum + x(i)   
       ENDDO
       xmean = sum / n 
    end subroutine calcmean  
    
    subroutine calcstd
!
! calculate the standard deviation
!                
      sum = 0.0
      DO i = 1 , n 
        sum = sum + (x(i)-xmean)**2
      ENDDO
      xstd = sqrt ( sum /(n-1 ) ) 
    end subroutine calcstd 
  
 
  END
