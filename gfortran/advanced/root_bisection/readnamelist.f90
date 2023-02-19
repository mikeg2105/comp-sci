	program readnamelist
	CHARACTER*14 SAMPLE 
	LOGICAL*4 NEW 
	REAL*4 DELTA, MAT(2,2) 
	NAMELIST /CASE/ SAMPLE, NEW, DELTA, MAT 
	READ ( 1, CASE )
	write(*,*) delta
	end
