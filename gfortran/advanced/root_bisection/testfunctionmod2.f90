module testfunctionmod
!  test funcction for a function of the form
contains


real function testfunction(x)


real, intent(in) ::  x


	testfunction=x**3-2*x-5


end function testfunction

end module testfunctionmod
