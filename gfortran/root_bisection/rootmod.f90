module rootmod


use testfunctionmod


implicit none

contains


real function bisection(a,b,fa,fb,eps,inf)

      real, intent(in)  :: eps
      real, intent(in) :: inf
      real, intent(inout)    ::  a , b , fa, fb
      real  :: x, test, fx

      test=1
      do while (abs(b-a)>(eps*b))

  	x = (a+b)/2
  	fx = testfunction(x)
  	if (sign(test,fx)==sign(test,fa)) then
    		a = x
    		fa = fx
    		print *, "left"
  	else
    		b = x
    		fb = fx
    		print *, "right"
  	endif
  	print * ,x,"  ",fx

      end do
      bisection=x




end function bisection


real function newton(a,b,fa,fb,eps,inf)

      real, intent(in)  :: eps
      real, intent(in) :: inf
      real, intent(inout)    ::  a , b , fa, fb
      real  :: x, test, fx, dfx

      test=1

      x=b
      do while ((sqrt((x-a)**2))>eps)
        a=x
        fx = testfunction(x)
        dfx= testderivfunction(x)
        x=a-(fx/dfx)
      
        print *,x-a
      end do


      newton=x


end function newton

! regula falsi method - linear interpolation method
real function interp(a,b,fa,fb,eps,inf)

      real, intent(in)  :: eps
      real, intent(in) :: inf
      real, intent(inout)    ::  a , b , fa, fb
      real  :: x, test, fx, dfx

      test=1
      x=b
      do while (abs(b-a)>(eps))


        fb=testfunction(b)
        fa=testfunction(a)
  	x=b-fb*(b-a)/(fb-fa)
  	fx = testfunction(x)
  	print *,a," ",fa," ",b," ",fb  	
  	
  	
  	if (sign(test,fx)==sign(test,fa)) then
    		b = x
    		fb = fx
    		print *, "left"
  	else
    		a = x
    		fa = fx
    		print *, "right"
  	endif
  	print * ,x,"  ",fx

      end do

      interp=x

end function interp

end
