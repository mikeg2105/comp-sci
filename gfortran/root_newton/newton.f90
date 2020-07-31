      program newton

! This is an example use of the plot function 
! used in association with the finding of a root
! problem solved earlier.
! % a is the lower value of the range
! b is the upper value of the range

      real, parameter :: eps = 1.0E-5
      real, parameter :: inf = 1.0E30
      real, parameter :: pi = 4*atan(1.0)
      real    ::  a , b , fa, fb, x, test, delta

      test=1
      a = 0
      fa = -inf
      x = 3
      fb = inf
      b=x

      delta = (b-a)/100

      do while ((sqrt((x-a)**2))>eps)
        a=x
        fx = 0.5*sin(2*(x-(pi/4)))+0.5*sin(x)
        dfx= cos(2*(x-(pi/4)))+0.5*cos(x)
        x=a-(fx/dfx)
      
        print *,x-a
      end do

      print *, ' The root is :',x
      

      end program 
