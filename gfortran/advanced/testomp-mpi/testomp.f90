program Parallel_Hello_World



! gfortran testomp.f90 -o testomp -fopenmp

! ./testomp
! echo $OMP_NUM_THREADS

USE OMP_LIB

implicit none

integer, parameter :: n = 4
integer :: i,ntot
real, dimension(n) :: dat, result


 ntot= OMP_GET_NUM_THREADS()

!$OMP PARALLEL
    dat=OMP_GET_THREAD_NUM()
    PRINT *, 'Hello from process: ', OMP_GET_THREAD_NUM(), "   ", i
    result=dat+1
!$OMP END PARALLEL

print *, result, " ",dat

END
