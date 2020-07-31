      program findrt


! here is an example use of the while statement
! which is used for finding the root of a polynomial 
! which is known to lie within a certain interval.
! a is the lower value of the range
! b is the upper value of the range
      real, parameter :: eps = 1.0E-5
      real, parameter :: inf = 1.0E30
      real    ::  a , b , fa, fb, x, test

      test=1
      a = 0
      fa = -inf
      b = 3
      fb = inf

      do while (abs(b-a)>(eps*b))

  	x = (a+b)/2
  	fx = x**3-2*x-5
  	if (sign(test,fx)==sign(test,fa)) then
    		a = x
    		fa = fx
  	else
    		b = x
    		fb = fx
  	endif

      end do

      print *, ' The root is :', x
      
      end
