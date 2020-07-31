! Exercise 5: Taylors series expansion of sin(x).
! 
!      sinx = x - x^3/3!  + x^5/5! - x7/7! + ......   
! terminate calculations when the next term to calculate becomes less than 1.E-5
! 
!
        program taylor     
        real, parameter :: epsilon = 1.0E-5
        integer ::  i 
        real    ::  angle=1 , pi , nexterm
                            
        print *, ' Enter an angle in degrees, enter 0 to terminate program:'
        
        do while ( angle /= 0 ) 
         read *, angle   
         
! convert angle to radians
         pi = 4*atan(1.0)  ! arc tangent of 1.0 is 45 degrees=1/4 radians.
         angle = angle*pi/180 

         sum = 0.0   
         nexterm = angle
         i = 1    
         do while ( abs(nexterm) >= epsilon) 
           sum = sum + nexterm
           nexterm = - nexterm*(angle**2)/( (i+1)*(i+2) )
           i = i + 2 
         end do
         
         sum = sum + nexterm          
         
         print *, sum , sin(angle) 
       enddo
       print * , 'finished...' 
        end
        
 