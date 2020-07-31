! This is another example solution for Exercise 1(b)
!
!
        program ex1b
        real :: a, b 
        real sum , prod , quo 
        print*,' Please enter two numbers'
        read *, a , b 
          sum = a + b ; prod = a*b ; quo = a/b 
          print *,'The sum of these numbers is:',sum
          print *,' Product is:' , prod
          print *,'Quotient is:' , quo 
        end 
