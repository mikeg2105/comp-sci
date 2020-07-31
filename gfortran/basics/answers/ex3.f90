!
! An example solution to exercise 3
! Calculates the area of a triangle with given length of three sides. 
!
       program areacalculate
       real  :: sidea , sideb , sidec , area , S , sum
       
         print*, 'Please Enter the length of the three sides of the triangle:'
         read *, sidea , sideb , sidec 
         
         sum =  sidea+sideb+sidec 
! draw a very acute triangle to visualise why the below condition is needed.  
         if( 2.0*max(sidea,sideb,sidec) >= sum ) then 
           print *, 'Invalid Triange !' 
         else 
           S = 0.5*sum 
           area = sqrt( S*(S-sidea)*(S-sideb)*(S-sidec) ) 
           print*,'Area of the triangle is' , area
         endif
         end
         