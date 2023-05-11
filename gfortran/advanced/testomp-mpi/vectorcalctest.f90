program simple

use omp_lib

! gfortran testomp.f90 -o testomp -fopenmp

! ./testomp
! echo $OMP_NUM_THREADS



!integer,parameter :: n=huge(1)
integer,parameter :: n=999999,n2=10000
real :: array(n)
integer :: iin(n)
real :: t1, t2, tt
print *,n

do i = 1, n
   iin(i) = i
enddo

print *,"loop  with mp - root an integer"
call cpu_time( t1 )
!$omp parallel do
do j=1,n2
do i = 1, n
   array(i) = sqrt(real(i))
enddo
enddo
!$omp end parallel do
call cpu_time( t2 )
tt=t2-t1
print *,tt



print *,"loop  with mp"
call cpu_time( t1 )
!$omp parallel do
do j=1,n2
do i = 1, n
   array(i) = sqrt(real(iin(i)))
enddo
enddo
!$omp end parallel do
call cpu_time( t2 )
tt=t2-t1
print *,tt

!print *,"array index  with mp"
!call cpu_time( t1 )
!!$omp parallel do
!   array(1:n) = sqrt(real(iin(1:n)))
!!$omp end parallel do
!call cpu_time( t2 )
!tt=t2-t1
!print *,tt

print *,"array index  no mp"
call cpu_time( t1 )
do j=1,n2
array(1:n) = sqrt(real(iin(1:n)))
enddo
call cpu_time( t2 )
tt=t2-t1
print *,tt

print *,"loop no mp"
call cpu_time( t1 )
do j=1,n2
do i = 1, n
   array(i) = sqrt(real(iin(i)))
enddo
enddo
call cpu_time( t2 )
tt=t2-t1
print *,tt






end program
