! Exercise 5: Taylors series expansion of sin(x).
! 
!      sinx = x - x^3/3!  + x^5/5! - x7/7! + ......  total of 7 terms.
!
        program taylor 
        integer ::  i 
        real    ::  angle , pi , nexterm
                            
        print *, ' Enter an angle in radians:'
        
         read *, angle   
         
! convert angle to radians
         pi = 4*atan(1.0)  ! arc tangent of 1.0 is 45 degrees=1/4 radians.
         angle = angle*pi/180 

         sum = 0.0   
         nexterm = angle    
         do i = 1 , 13 , 2  
           sum = sum + nexterm
           nexterm = - nexterm*(angle**2)/( (i+1)*(i+2) ) 
         enddo          
         
         print *, sum , sin(angle) 
        
        end
        
 