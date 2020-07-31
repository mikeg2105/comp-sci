  PROGRAM stats
!
! program to calculate the mean and standard deviation of a set of data values.
! 
  IMPLICIT NONE
  INTEGER, PARAMETER :: NMAX=200 
  INTEGER :: n           
  REAL :: x(NMAX) , xmean , xstd 
  INTERFACE
    subroutine readdata( a, n)
       real , intent(out) :: a(:)
       integer , intent(out) :: n
    end subroutine readdata
    real function calcmean( a , n )                                     
      real, intent(in) :: a(:) 
      integer, intent(in) :: n 
    end function calcmean
    real function calcstd  ( a ,xmean,  n ) 
      real, intent(in) :: a(:) , xmean 
      integer, intent(in) :: n 
    end
  END INTERFACE
  
  call readdata( x, n )
  xmean = calcmean( x, n )
  xstd = calcstd( x ,xmean, n )
 
  write (*,100) xmean , xstd 
100 format ( 'MEAN VALUE=' , F12.4 ,5X,'STANDARD DEVIATION=',F12.6 )
  stop
  end
  
  subroutine readdata(a , n )
  IMPLICIT NONE
  real, intent(out) :: a(:) 
  integer , intent(out) :: n 
  integer :: i 
  
    print * , 'Enter the number of data values to be read'  
    read *, n
    print *, 'Enter the data values:' 
    read (*,*) ( a(i) , i = 1 , n )  
  end  
    
  real function calcmean( a , n )                               
  IMPLICIT NONE        
  real, intent(in) :: a(:) 
  integer, intent(in) :: n 
  integer :: i  
  real   :: sum 
!
! calculate the mean
!     
       sum = 0.0              
       DO i = 1  , n
         sum = sum + a(i)   
       ENDDO
       calcmean = sum / n 
   end 
    
   real function calcstd  ( a ,xmean,  n ) 
   IMPLICIT NONE   
   real, intent(in) :: a(:) , xmean 
   integer, intent(in) :: n 
   integer :: i 
   real   :: sum     
!
! calculate the standard deviation
!                
      sum = 0.0
      DO i = 1 , n 
        sum = sum + (a(i)-xmean)**2
      ENDDO
      calcstd =sqrt(  sum /(n-1 ) ) 
   end  
  
