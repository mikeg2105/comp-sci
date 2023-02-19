      program findrt

      use testfunctionmod
      use rootmod
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
      b = 5
      fb = inf
      
      !use the bisection function in the rootmod module
      x=interp(a,b,fa,fb,eps,inf)
      
      print *, ' The root is :', x


      
      end
