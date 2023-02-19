module testfunctionmod
!  test funcction for a function of the form



contains


real function testfunction(x)


real, intent(in) ::  x


	testfunction=x**3+x**2-3*x-3


end function testfunction

real function testderivfunction(x)


real, intent(in) ::  x


	testderivfunction=3*x**2+2*x-3


end function testderivfunction


end module testfunctionmod
