  PROGRAM stats
!
! program to calculate basic statistics of data in a file.
! 
  IMPLICIT NONE
  INTEGER, PARAMETER :: NMAX=200 
  INTEGER :: n , i    
  character*64 :: fname      
  REAL :: x(NMAX) , xmean , xstd , sum
  
10  print * , 'Enter the name of your datafile:'  
    read(*,120) fname
120 format( A ) 
    OPEN ( UNIT = 8, FILE=fname ,STATUS='OLD',ERR=10 )
  
    do i = 1 , NMAX
      read(8,*, end=30) x(i)
    enddo
    i= NMAX+1
30  n =i-1
    print *,'Number of data points=' , n
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
! open an output file    
   open ( unit=9, FILE='normalised.txt' , STATUS='UNKNOWN' )
! write normalised values into it. 
   write(9,102) 
102 format( ' Normalised Values:' ) 
   DO i = 1 , n 
      write (9,104)( x(i)-xmean ) / xstd   
   ENDDO
104 format (G18.6 )  
 
   close(8)
   close (9) 
   
  END
