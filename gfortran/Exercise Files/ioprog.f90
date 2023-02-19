program io
implicit none
character(len=20) :: filename
real :: x, y
filename='results.out'
x = 1.
y = 2.
! write to ascii file
open(unit=1,file=filename,status='replace',form='formatted')
write(1,*) x, y
close(1)
! reset vars
x = 0.
y = 0.
! read from ascii file
open(unit=2,file=filename,status='old')
read(2,*) x,y
close(2)
! print vars
print*,' x = ',x,' y = ',y
end program io
