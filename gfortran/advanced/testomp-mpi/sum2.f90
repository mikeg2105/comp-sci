      PROGRAM sum2
!
!     Program to sum all elements of a two dimensional array.
!     Version 2: Fortran 77
!                Correct elements passed by defining array in
!                subroutine to have first dimension the same as
!                in the main program.

! gfortran sum2.f90 -o testomp -fopenmp

! ./testomp
! echo $OMP_NUM_THREADS
USE OMP_LIB

      PARAMETER (n=10)
      REAL ra(n,n)
      INTEGER i, j, n1, nsum

    integer :: ntot
    real, dimension(n) :: dat, result

    ntot= OMP_GET_NUM_THREADS()

    !$OMP PARALLEL
        dat=OMP_GET_THREAD_NUM()
        PRINT *, 'Hello from process: ', OMP_GET_THREAD_NUM(), "   ", i
        result=dat+1
    !$OMP END PARALLEL



    !$OMP DO
      DO 20 i=1,n
        DO 10 j=1,n
          ra(j,i) =1000
   10   CONTINUE
   20 CONTINUE
   !$OMP ENDDO

      WRITE (*,'("Input matrix size:")')
      READ (*,*) n1

    !$OMP DO
      DO 40 i=1,n1
        DO 30 j=1,n1
          ra(j,i) = (i-1)*n1 + j
   30   CONTINUE
   40 CONTINUE
   !$OMP ENDDO

      CALL sumsub(ra,n1,nsum,n)
      WRITE (*,'("Sum for matrix of size ",I3," = ",I10)') n1,nsum

      END

      SUBROUTINE sumsub(ra,n1,nsum,n)

      REAL ra(n,n1)
      INTEGER n1, nsum, i, j

      nsum = 0
      !$OMP DO
      DO 20 i=1,n1
        DO 10 j=1,n1
          nsum = nsum + ra(j,i)
   10   CONTINUE
   20 CONTINUE
   !$OMP ENDDO

      RETURN
      END
