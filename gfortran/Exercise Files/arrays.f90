module arrays
implicit none
contains
subroutine array_examples()
use prec, only:sp
real(sp) :: x(3),y(3)
! set each element individually
x(1) = 1.0
x(2) = 2.0
x(3) = 3.0
print*,' x = ',x
! set whole array equal to 1
x = 1.
print*,' x = ',x
! set array parts in one line
x = (/1.0, 2.0, 3.0/)
print*,' x = ',x
! set y=x
y = x
print*,' y = ',y
! array operations
print*,'sum=',sum(y)
print*,'maxval=',maxval(y)
print*,'minval=',minval(y)
print*,'maxloc=',maxloc(y)
print*,'minloc=',minloc(y)
! magnitude
print*,' |x| = ',sqrt(dot_product(x,x))
print*,' |x| = ',norm2(x)
end subroutine array_examples
end module arrays
