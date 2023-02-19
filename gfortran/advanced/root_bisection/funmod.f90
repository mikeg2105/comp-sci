module testfunctionmod
!  test funcction for a function of the form
contains


real function testfunction(x)


real, intent(in) ::  x
real, intent(out) :: y

	y=x**3-2*x-5


end function testfunction

end module testfunctionmod
