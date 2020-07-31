      program interpol
!Integration using Simpsons rule

!quadrature to higher order using lagrange interpolation

      real, parameter :: eps = 1.0E-5
      real, parameter :: inf = 1.0E30
      real, parameter :: pi = 4*atan(1.0)
      integer, parameter :: n = 10
      integer, parameter :: nmax = (n*2)+1
      real    ::  a , b ,h, x(nmax),y(nmax), test
      integer :: i






!integration interval
      a=0
      b=pi/2

      h=(b-a)/n

     

      do i=1,nmax+1
         x(i)=a+h*(i-1)
         y(i)=testfunction(x(i))
      end do

      sum=0

!since we are using interpolation
!the loop performs the summation 
!over every other element
      do i=2,n-1,2
         sum=sum+y(i-1)+4*y(i)+y(i+1)
      end do

      result=sum*h/3

      print *,'results is',result

      end program


      real function testfunction (x)

      IMPLICIT NONE        
      real, intent(in) :: x 
 
      testfunction=sin(x)
      end

!lagrange interpolation

!xval is the value of x at which we wish to interpolate
!f is the series of values
      real function lagrangeinterp(xval,x,f,i)
      IMPLICIT NONE        
      real, intent(in) :: x(:),f(:),xval 
      integer, intent(in) :: i 
      real :: t1,t2,t3
                
      t1=(xval-x(i))*(xval-x(i+1))/((x(i-1)-x(i))*(x(i-1)-x(i+1)))
      t2=(xval-x(i-1))*(xval-x(i+1))/((x(i)-x(i-1))*(x(i)-x(i+1)))
      t3=(xval-x(i-1))*(xval-x(i))/((x(i+1)-x(i-1))*(x(i+1)-x(i)))
      lagrangeinterp=t1*f(i-1)+t2*f(i)+t3*f(i+1)
      end

