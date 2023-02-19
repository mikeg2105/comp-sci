      PROGRAM SHALLOW
      implicit none
C ***
C *** BENCHMARK WEATHER PREDICTION PROGRAM FOR COMPARING THE PERFORMANCE
C *** OF CURRENT SUPERCOMPUTERS. THE MODEL IS BASED ON THE PAPER, "THE
C *** DYNAMICS OF FINITE-DIFFERENCE MODELS OF THE SHALLOW-WATER
C *** EQUATIONS," BY R. SADOURNY, J. ATM. SCI., VOL.32, NO.4, APRIL 1975
C *** CODE BY PAUL N. SWARZTRAUBER, NATIONAL CENTER FOR ATMOSPHERIC
C *** RESEARCH, BOULDER, COLORADO   OCTOBER 1984.
*
* Converted to F90 array syntax - CA 31/07/2003
* Modified to use MPI           - CA 01/08/2003
* Checkpoint / restart added    - CA 11/08/2003
*
C ***

      integer n,m,itmax, mp1, np1, mnmin, ncycle, cycle_start
      integer i,j, mprint
      integer problem_size
      double precision   dx, dy, dt,tdt,a,alpha, time, ptime
      double precision,allocatable :: P(:,:),U(:,:),V(:,:),PSI(:,:)
      double precision,allocatable :: Pold(:,:),Uold(:,:),Vold(:,:)
      double precision,allocatable :: Pnew(:,:),Unew(:,:),Vnew(:,:)
      double precision,allocatable :: CU(:,:),CV(:,:),Z(:,:),H(:,:)
      double precision pi, di, dj, tpi, fsdx, fsdy,
     1                  tdts8, tdtsdx, tdtsdy, pcf
      double precision secs, clock, dsecnd
      double precision TCYC, XFLOPS
      double precision TOTTIM,CTIME
      double precision PMAX,PMIN,PAVE,PMAX_local,PMIN_local,PAVE_local
      integer imax, jmax, imin, jmin
      integer flush, system
      integer checkpoint, restart
      include "mpif.h"
      integer nx, ny
      integer myid, numprocs, it, rc, comm2d, ierr, stridetype
      integer nbrleft, nbrright, nbrtop, nbrbottom
      integer topLeft, topRight, bottomLeft, bottomRight
      integer sx, ex, sy, ey
      integer dims(2)
      logical, parameter:: periods(2) = (/.true., .true./)
      logical:: diag = .false.

      double precision diff2d, diffnorm, dwork
      double precision t1, t2
      integer status(MPI_STATUS_SIZE)
      external diff2d
      character(len=*),parameter:: root='./'
      character(len=80):: directory= root // 'data_store'
      character(len=80):: base_filename='savefile'
      character(len=*),parameter:: input_file='shall_input'
      character(len=4) myid_str, nproc_str, int_to_string
      integer,parameter:: chkpnt_unit=16
      integer,parameter:: restart_unit=15
      double precision, parameter:: base_pressure=5.d4
      integer buffer(6)     
      data checkpoint/0/
      data restart/0/

      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      print *, "Process ", myid, " of ", numprocs, " is alive"
*     open(5,file=input_file)
*
* These probably should be set-up as parameters
*     

      DT = 90.D0
      TDT = 90.D0
      TIME = 0.D0
      DX = 1.D5
      DY = 1.D5
      ALPHA = .001D0


      if (myid .eq. 0) then
       write(*,*) 'Input the problem size and number of time steps '
*       read(5,*) m,n,itmax,mprint,restart,checkpoint
       read(5,*) (buffer(i),i=1,6)
      endif
*
* Static restart information should be read here.
*     
      call mpi_bcast(buffer,6,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

*
*     m and n are global sizes,
*     nx and ny are the corresponding local sizes
*     
      m = buffer(1)
      n = buffer(2)
      itmax = buffer(3)
      mprint = buffer(4)
*
* restart at step x+1 assuming data read in is for step x if restart >  0
*   
      restart = buffer(5)
      checkpoint = buffer(6)
      if (myid .eq. 0) then
C *** PRINT INITIAL VALUES
          WRITE(*,390) N,M,DX,DY,DT,itmax,mprint,restart,checkpoint
  390     FORMAT(' NUMBER OF POINTS IN THE X DIRECTION',I8/
     1       ' NUMBER OF POINTS IN THE Y DIRECTION',I8/
     2       ' GRID SPACING IN THE X DIRECTION    ',F8.0/
     3       ' GRID SPACING IN THE Y DIRECTION    ',F8.0/
     4       ' TIME STEP                          ',F8.0/
     5       ' Max time steps                     ',i5/
     5       ' Check progress frequency           ',i5/
     6       ' Restart                            ',i5/
     7       ' Checkpoint at                      ',i5)

       end IF   
c
c Get a new communicator for a decomposition of the domain.  Let MPI
c find a "good" decomposition but all processes are used and the processes
c are not renumbered!
c
      if (restart .le. 0) then
         dims(1) = 0
         dims(2) = 0
         call MPI_DIMS_CREATE( numprocs, 2, dims, ierr )
      else
         open(restart_unit,form='unformatted',
     1     file=trim(directory)//'/'//trim(base_filename)//
     2          int_to_string(myid)//'of'//int_to_string(numprocs))
         read(restart_unit) dims
      end if
      call MPI_CART_CREATE( MPI_COMM_WORLD, 2, dims, periods, .false.,
     1                    comm2d, ierr )

c      print *, "Process ", myid, " of ", numprocs, " is alive"
c
c My neighbors are now +/- 1 with my rank. 
c
      call fnd2dnbrs( comm2d, nbrleft, nbrright, nbrtop, nbrbottom )
      topLeft = 0
      bottomRight = dims(1)*dims(2)-1
      topRight = dims(2) - 1
      bottomLeft = bottomRight - dims(2) + 1
      if (diag) then
        print *, "Process ", myid, " of ", numprocs, " is alive"
        print *, "Process ", myid, ":", 
     *     nbrleft, nbrright, nbrtop, nbrbottom
      end if
      if (myid .eq. 0 .and. diag) then
        print *, 'Corners ',topLeft, topRight, bottomLeft, bottomRight
      end if 
c
c Compute the decomposition
c     
      call fnd2ddecomp( comm2d, m, n, sx, ex, sy, ey )
      print *, "Process ", myid, ":", sx, ex, sy, ey
      nx = ex - sx + 1
      ny = ey - sy + 1
c
c Create a new, "strided" datatype for the exchange in the "non-contiguous"
c direction
c
      call mpi_Type_vector( ey-sy+1, 1, ex-sx+3, 
     $                      MPI_DOUBLE_PRECISION, stridetype, ierr )
      call mpi_Type_commit( stridetype, ierr )       

      allocate (P(sx-1:ex+1,sy-1:ey+1))
      allocate (U(sx-1:ex+1,sy-1:ey+1))
      allocate (V(sx-1:ex+1,sy-1:ey+1))
      allocate (Psi(sx-1:ex+1,sy-1:ey+1))   
      allocate (Pold(sx-1:ex+1,sy-1:ey+1))
      allocate (Uold(sx-1:ex+1,sy-1:ey+1))
      allocate (Vold(sx-1:ex+1,sy-1:ey+1))   
      allocate (Pnew(sx-1:ex+1,sy-1:ey+1))
      allocate (Unew(sx-1:ex+1,sy-1:ey+1))
      allocate (Vnew(sx-1:ex+1,sy-1:ey+1))   
      allocate (CU(sx-1:ex+1,sy-1:ey+1))   
      allocate (CV(sx-1:ex+1,sy-1:ey+1))   
      allocate (Z(sx-1:ex+1,sy-1:ey+1))   
      allocate (H(sx-1:ex+1,sy-1:ey+1))
      if (restart .le. 0) then
        call init()
        cycle_start = 1

      else
        read(restart_unit) cycle_start
        read(restart_unit) p(sx:ex+1,sy:ey+1)
        read(restart_unit) pold(sx:ex+1,sy:ey+1)
        read(restart_unit) u(sx:ex+1,sy:ey+1)
        read(restart_unit) uold(sx:ex+1,sy:ey+1)
        read(restart_unit) v(sx:ex+1,sy:ey+1)
        read(restart_unit) vold(sx:ex+1,sy:ey+1)
        close(restart_unit)
        time = (cycle_start)*DT
*
* The above data values represent the completed step cycle_start
* Therefore, the next step should be one higher. Also as at least one
* step has been completed, TDT should be double the value of DT.
*       
        cycle_start = cycle_start + 1
        TDT = DT + DT
*       
      end if



!      clock = dsecnd()

      do ncycle = cycle_start, itmax

         call capuvhz()

C *** COMPUTE NEW VALUES U, V, AND P
      TDTS8 = TDT/8.d0
      TDTSDX = TDT/DX
      TDTSDY = TDT/DY


      DO J=sy,ey
         DO I=sx,ex
            UNEW(I+1,J) = UOLD(I+1,J)+
     1        TDTS8*(Z(I+1,J+1)+Z(I+1,J))*(CV(I+1,J+1)+CV(I,J+1)+CV(I,J)
     2       +CV(I+1,J))-TDTSDX*(H(I+1,J)-H(I,J))
            VNEW(I,J+1) = VOLD(I,J+1)-TDTS8*(Z(I+1,J+1)+Z(I,J+1))
     1       *(CU(I+1,J+1)+CU(I,J+1)+CU(I,J)+CU(I+1,J))
     2       -TDTSDY*(H(I,J+1)-H(I,J))
            PNEW(I,J) = POLD(I,J)-TDTSDX*(CU(I+1,J)-CU(I,J))
     1       -TDTSDY*(CV(I,J+1)-CV(I,J))
         END DO
      END DO
      call MPI_SENDRECV(unew(ex+1,sy),1,stridetype,nbrbottom, 0,
     1                  unew(sx,sy),1,stridetype,nbrtop, 0,
     2                  comm2d, status, ierr)
      call MPI_SENDRECV(vnew(sx,sy+1),1,stridetype,nbrtop, 1,
     1                  vnew(ex+1,sy+1),1,stridetype,nbrbottom, 1,
     2                  comm2d, status, ierr)
      call MPI_SENDRECV(pnew(sx,sy),1,stridetype,nbrtop, 2,
     1                  pnew(ex+1,sy),1,stridetype,nbrbottom, 2,
     2                  comm2d, status, ierr)
C *** PERIODIC CONTINUATION - Needs communication
*      DO 210 J=1,N
*      UNEW(1,J) = UNEW(M+1,J)
*      VNEW(M+1,J+1) = VNEW(1,J+1)
*      PNEW(M+1,J) = PNEW(1,J)
*  210 CONTINUE
      call MPI_SENDRECV(unew(sx+1,sy),nx,MPI_DOUBLE_PRECISION,nbrleft,3,
     1                 unew(sx+1,ey+1),nx,MPI_DOUBLE_PRECISION,
     2                 nbrright,3, comm2d, status, ierr)
      call MPI_SENDRECV(vnew(sx,ey+1),nx,MPI_DOUBLE_PRECISION,
     1                  nbrright,4,
     2                  vnew(sx,sy),nx,MPI_DOUBLE_PRECISION,
     3                  nbrleft,4,comm2d, status, ierr)
      call MPI_SENDRECV(pnew(sx,sy),nx,MPI_DOUBLE_PRECISION,nbrleft, 5,
     1                 pnew(sx,ey+1),nx,MPI_DOUBLE_PRECISION,nbrright,5,
     2                 comm2d, status, ierr)
*      DO 215 I=1,M
*      UNEW(I+1,N+1) = UNEW(I+1,1)
*      VNEW(I,1) = VNEW(I,N+1)
*      PNEW(I,N+1) = PNEW(I,1)
*  215 CONTINUE
* Top Right     
      if (myid .eq. topRight .and. myid .ne. bottomLeft) then
        call MPI_RECV(unew(sx,ey+1),1, MPI_DOUBLE_PRECISION, bottomLeft,
     1                 m,comm2d,status,ierr)
        call MPI_SEND(vnew(sx,ey+1),1, MPI_DOUBLE_PRECISION,bottomLeft,
     1                 m+1,comm2d,ierr)
      end if
* Bottom Left     
      if (myid .eq. bottomLeft .and. myid .ne. topRight) then
        call MPI_SEND(unew(ex+1,sy),1, MPI_DOUBLE_PRECISION,topRight,
     1                 m,comm2d,ierr)
        call MPI_RECV(vnew(ex+1,sy),1, MPI_DOUBLE_PRECISION, topRight,
     1                 m+1,comm2d,status,ierr)
      end if
      if (myid .eq. bottomLeft .and. myid .eq. topRight) then
       UNEW(1,N+1) = UNEW(M+1,1)
       VNEW(M+1,1) = VNEW(1,N+1) 
      end if   
*Top Left     
      if (myid .eq. topLeft .and. myid .ne. bottomRight) then

         call MPI_SEND(pnew(sx,sy),1, MPI_DOUBLE_PRECISION,bottomRight,
     1                 m,comm2d,ierr)
      end if
*Bottom Right     
      if (myid .eq. bottomRight .and. myid .ne. topLeft) then

         call MPI_RECV(pnew(ex+1,ey+1),1, MPI_DOUBLE_PRECISION, topLeft,
     1                 m,comm2d,status,ierr)
      end if
      if (myid .eq. bottomRight .and. myid .eq. topLeft) then
       PNEW(M+1,N+1) = PNEW(1,1)
      end if

*
* Note - communications are NOT required here because the copied elements
* ( u, uNEW etc. have already been copied to match periodic bdy conds.)
*
      IF(NCYCLE .gt. 1) then
       DO J=sy,ey+1
         DO I=sx,ex+1
            UOLD(I,J) = U(I,J)+ALPHA*(UNEW(I,J)-2.*U(I,J)+UOLD(I,J))
            VOLD(I,J) = V(I,J)+ALPHA*(VNEW(I,J)-2.*V(I,J)+VOLD(I,J))
            POLD(I,J) = P(I,J)+ALPHA*(PNEW(I,J)-2.*P(I,J)+POLD(I,J))
            U(I,J) = UNEW(I,J)
            V(I,J) = VNEW(I,J)
            P(I,J) = PNEW(I,J)
         END DO
       END DO
      else
       TDT = TDT+TDT
       DO J=sy,ey+1
         DO I=sx,ex+1
            UOLD(I,J) = U(I,J)
            VOLD(I,J) = V(I,J)
            POLD(I,J) = P(I,J)
            U(I,J) = UNEW(I,J)
            V(I,J) = VNEW(I,J)
            P(I,J) = PNEW(I,J)
         END DO
       END DO
      end if
      TIME = TIME + DT
*
* Checkpointing should take place here and should involve these 6 arrays
* plus some static problem definition info written initially.
*       
      if (checkpoint .gt. 0 .and. mod(ncycle,checkpoint) .eq. 0) then
         ierr = system('mkdir  -p '//trim(directory))
         open(chkpnt_unit,form='unformatted',
     1          file=trim(directory)// '/'//trim(base_filename)//
     2          int_to_string(myid)//'of'//int_to_string(numprocs))
         write( chkpnt_unit) dims
         write(chkpnt_unit) ncycle
         write(chkpnt_unit) p(sx:ex+1,sy:ey+1)
         write(chkpnt_unit) pold(sx:ex+1,sy:ey+1)
         write(chkpnt_unit) u(sx:ex+1,sy:ey+1)
         write(chkpnt_unit) uold(sx:ex+1,sy:ey+1)
         write(chkpnt_unit) v(sx:ex+1,sy:ey+1)
         write(chkpnt_unit) vold(sx:ex+1,sy:ey+1)   
         close(chkpnt_unit)
      end if
      IF(MOD(NCYCLE,MPRINT).eq.0 .or. NCYCLE .eq. ITMAX) then
         PTIME = TIME/3600.

*         WRITE(*,355) (PNEW(I,I),I=1,MNMIN)
*  355    FORMAT(/,' DIAGONAL ELEMENTS OF P ' //,(8E16.6))
*         WRITE(*,360) (UNEW(I,I),I=1,MNMIN)
*  360    FORMAT(/,' DIAGONAL ELEMENTS OF U ' //,(8E16.6))
*         WRITE(*,365) (VNEW(I,I),I=1,MNMIN)
*  365    FORMAT(/,' DIAGONAL ELEMENTS OF V ' //,(8E16.6))

C        TIME IN SECONDS
!         tottim = max(1.d-6,(dsecnd() - clock))

C        COMPUTE THE MFLOPS FOR THE RUN
!         XFLOPS = float(N*M)*float(NCYCLE)*65*1.0e-6/TOTTIM

C-----------------------------------------------------------
C        COMPUTE MIN,MAX AND SUM VALUES FOR DIAGNOSTICS
C-----------------------------------------------------------
         PMIN=1.D30
         PMAX=-1D30
         PAVE=0.0D0

         do j=sy,ey
            do i=sx,ex
               PMAX=MAX(PMAX,P(i,j))
*               if (p(i,j) .gt. pmax) then
*                 pmax = p(i,j)
*                 imax = i
*                 jmax = j
*              end if
               PMIN=MIN(PMIN,P(i,j))
*               if (p(i,j) .lt. pmin) then
*                 pmin = p(i,j)
*                 imin = i
*                 jmin = j
*              end if
               PAVE=PAVE+P(i,j)
            end do
         end do
*        PMAX=PMAX-50000.
*        PMIN=PMIN-50000.
*        PAVE=PAVE/(M*N)-50000.
         PMAX_local = PMAX
         PMIN_local = PMIN
         PAVE_local = PAVE
         secs = tottim
         if (diag) then
           write(*,370) ncycle,(j,p(j,j)-base_pressure,
     1                          j=max(sx,sy),min(ex,ey))
  370      format(i4,500(i4,1pd22.15)) 
         end if
         call mpi_reduce(PMAX_local,PMAX,1,MPI_DOUBLE_PRECISION,
     1                   MPI_MAX,0,comm2d,ierr)
         call mpi_reduce(PMIN_local,PMIN,1,MPI_DOUBLE_PRECISION,
     1                   MPI_MIN,0,comm2d,ierr)         
         call mpi_reduce(PAVE_local,PAVE,1,MPI_DOUBLE_PRECISION,
     1                   MPI_SUM,0,comm2d,ierr)         
         call mpi_reduce(secs,tottim,1,MPI_DOUBLE_PRECISION,
     1                   MPI_MAX,0,comm2d,ierr)
*          print *, myid,ncycle,' Local Maxloc ',imax,jmax
*          WRITE(6,9060)PMAX_local - 50000.
*          WRITE(6,9070)PMIN_local - 50000.
*          WRITE(6,9080)PAVE_local/(nx*ny) - 50000.
*          print *, myid,ncycle,' Local Minloc ',imin,jmin                               
*
* Need reduction functions to get values across nodes and print results
*
         if (myid .eq. 0) then
          WRITE(*,350) NCYCLE,PTIME
  350     FORMAT(//,' CYCLE NUMBER',I5,' MODEL TIME IN  HOURS',F6.2)
C        COMPUTE THE MFLOPS FOR THE RUN
          XFLOPS = float(N*M)*float(NCYCLE-cycle_start+1)*
     1            65*1.0e-6/TOTTIM
          if (diag) then
           do i=sx,sx+4
             write(*,9000) i,
     1                     (p(i,j)-base_pressure,j=sx,sx+5)
 9000        format(' P ',i3,6(2x,1pd10.3))
           end do
           do i=sx,sx+4
             write(*,9001) i,
     1                     (cu(i,j),j=sx,sx+5)
 9001        format('CU ',i3,6(2x,1pd10.3))
           end do       
           do i=sx,sx+4
             write(*,9002) i,
     1                     (cv(i,j),j=sx,sx+5)
 9002        format('CV ',i3,6(2x,1pd10.3))
           end do
          end if       
          PMAX=PMAX-base_pressure
          PMIN=PMIN-base_pressure
          PAVE=PAVE/(M*N)-base_pressure
          WRITE(6,9030) M, N
          print *, '-------------------------------------------------'
          WRITE(6,9060)PMAX
          WRITE(6,9070)PMIN
          WRITE(6,9080)PAVE
          print *, '-------------------------------------------------'
          WRITE(6,9040)TOTTIM
          WRITE(6,9050)XFLOPS
9030      FORMAT(' RESOLUTION = ',I5,' BY ',I5)
9040      FORMAT(' WALL CLOCK TIME FOR JOB = ', F10.5, ' seconds')
9050      FORMAT(' EXPECTED MFLOPS RATE    = ', F11.5)
*
*         Restrict decimal range to account for significant figures
*         This needs changing if the underlying problem is changed!!
*
9060      FORMAT(' MAX P FIELD             = ', E22.11)
9070      FORMAT(' MIN P FIELD             = ', E22.11)
9080      FORMAT(' AVG P FIELD             = ', E22.11)
         end if

      END IF     
      end do
C
      WRITE (*,220)
  220 FORMAT('0   *****  END OF PROGRAM SHALLOW  *****')
C

!      secs = dsecnd() - clock
!      write(*,*) 'Time taken: ',secs
        deallocate (P)
        deallocate (U)
        deallocate (V)
        deallocate (Psi) 
        deallocate (Pold)
        deallocate (Uold)
        deallocate (Vold)   
        deallocate (Pnew)
        deallocate (Unew)
        deallocate (Vnew)   
        deallocate (CU)   
        deallocate (CV)   
        deallocate (Z)   
        deallocate (H)
      call MPI_Type_free( stridetype, ierr )
      call MPI_Comm_free( comm2d, ierr )
      call MPI_FINALIZE(ierr)

      contains

      subroutine init()
      A = 1.D6         
      PI = 4.D0*ATAN(1.D0)
      TPI = PI+PI
      DI = TPI/dble(M)
      DJ = TPI/dble(N) 
*
* Init should not be executed on a restart
*     
      pcf=(pi**2)*(a**2)/((n*dx)**2)
C *** INITIAL VALUES OF THE STREAMFUNCTION

      DO J=sy,ey+1
         DO I=sx,ex+1
           PSI(I,J) = A*SIN((dble(I)-.5d0)*DI)*SIN((dble(J)-.5d0)*DJ)
           P(i,j)=pcf*(cos(2.d0*dble(i-1)*di)+
     1                 cos(2.d0*dble(j-1)*dj))+base_pressure
           POLD(I,J) = p(i,j)
         END DO
      END DO
C *** INITIALIZE VELOCITIES
      DO J=sy,ey
        DO I=sx,ex
           U(I+1,J) = -(PSI(I+1,J+1)-PSI(I+1,J))/DY
           V(I,J+1) = (PSI(I+1,J+1)-PSI(I,J+1))/DX
        end do
      end do
C *** PERIODIC CONTINUATION - Needs communication - Encapsulate later
*      DO 70 J=1,N
*      U(1,J) = U(M+1,J)
*      V(M+1,J+1) = V(1,J+1)
*   70 CONTINUE
      call MPI_SENDRECV(u(ex+1,sy),1,stridetype,nbrbottom, 0,
     1                  u(sx,sy),1,stridetype,nbrtop, 0,
     2                  comm2d, status, ierr)
      call MPI_SENDRECV(v(sx,sy+1),1,stridetype,nbrtop, 1,
     1                  v(ex+1,sy+1),1,stridetype,nbrbottom, 1,
     2                  comm2d, status, ierr)
*      DO 75 I=1,M
*      U(I+1,N+1) = U(I+1,1)
*      V(I,1) = V(I,N+1)
*   75 CONTINUE
      call MPI_SENDRECV(u(sx+1,sy),nx,MPI_DOUBLE_PRECISION,nbrleft, 0,
     1                 u(sx+1,ey+1),nx,MPI_DOUBLE_PRECISION,nbrright,0,
     2                 comm2d, status, ierr)
      call MPI_SENDRECV(v(sx,ey+1),nx,MPI_DOUBLE_PRECISION,nbrright,1,
     1                 v(sx,sy),nx,MPI_DOUBLE_PRECISION,nbrleft,1,
     2                 comm2d, status, ierr)
* Top Right - watch out for deadlock     
      if (myid .eq. topRight .and. myid .ne. bottomLeft) then
        call MPI_RECV(u(sx,ey+1),1, MPI_DOUBLE_PRECISION, bottomLeft,
     1                 m,comm2d,status,ierr)
        call MPI_SEND(v(sx,ey+1),1, MPI_DOUBLE_PRECISION,bottomLeft,
     1                 m+1,comm2d,ierr)
      end if
* Bottom Left     
      if (myid .eq. bottomLeft .and. myid .ne. topRight) then
        call MPI_SEND(u(ex+1,sy),1, MPI_DOUBLE_PRECISION, topRight,
     1                 m,comm2d, ierr)
        call MPI_RECV(v(ex+1,sy),1, MPI_DOUBLE_PRECISION, topRight,
     1                 m+1,comm2d,status,ierr)
      end if
* Single processor case     
      if (myid .eq. bottomLeft .and. myid .eq. topRight) then 
       U(1,N+1) = U(M+1,1)
       V(M+1,1) = V(1,N+1)
      end if
      DO J=sy,ey+1
         DO I=sx,ex+1
            UOLD(I,J) = U(I,J)
            VOLD(I,J) = V(I,J)
         END DO
      END DO
      end subroutine init

      subroutine capuvhz()
C *** COMPUTE CAPITAL U, CAPITAL V, Z, AND H
      FSDX = 4.D0/DX
      FSDY = 4.D0/DY
      DO J=sy,ey
         DO I=sx,ex
            CU(I+1,J) = .5d0*(P(I+1,J)+P(I,J))*U(I+1,J)
            CV(I,J+1) = .5d0*(P(I,J+1)+P(I,J))*V(I,J+1)
            Z(I+1,J+1) = (FSDX*(V(I+1,J+1)-V(I,J+1))-FSDY*(U(I+1,J+1)
     1          -U(I+1,J)))/(P(I,J)+P(I+1,J)+P(I+1,J+1)+P(I,J+1))
            H(I,J) = P(I,J)+.25d0*(U(I+1,J)*U(I+1,J)+U(I,J)*U(I,J)
     1               +V(I,J+1)*V(I,J+1)+V(I,J)*V(I,J))
         END DO
      END DO
C *** PERIODIC CONTINUATION - Needs communication
*      DO 110 J=1,N
*      CU(1,J) = CU(M+1,J)
*      CV(M+1,J+1) = CV(1,J+1)
*      Z(1,J+1) = Z(M+1,J+1)
*      H(M+1,J) = H(1,J)
*  110 CONTINUE
      call MPI_SENDRECV(cu(ex+1,sy),1,stridetype,nbrbottom, 0,
     1                  cu(sx,sy),1,stridetype,nbrtop, 0,
     2                  comm2d, status, ierr)
      call MPI_SENDRECV(cv(sx,sy+1),1,stridetype,nbrtop, 1,
     1                  cv(ex+1,sy+1),1,stridetype,nbrbottom, 1,
     2                  comm2d, status, ierr)
      call MPI_SENDRECV(z(ex+1,sy+1),1,stridetype,nbrbottom, 2,
     1                  z(sx,sy+1),1,stridetype,nbrtop, 2,
     2                  comm2d, status, ierr)
      call MPI_SENDRECV(h(sx,sy),1,stridetype,nbrtop, 3,
     1                  h(ex+1,sy),1,stridetype,nbrbottom, 3,
     2                  comm2d, status, ierr)     
*      DO 115 I=1,M
*      CU(I+1,N+1) = CU(I+1,1)
*      CV(I,1) = CV(I,N+1)
*      Z(I+1,1) = Z(I+1,N+1)
*      H(I,N+1) = H(I,1)
*  115 CONTINUE
      call MPI_SENDRECV(cu(sx+1,sy),nx,MPI_DOUBLE_PRECISION,nbrleft, 4,
     1                 cu(sx+1,ey+1),nx,MPI_DOUBLE_PRECISION,nbrright,4,
     2                 comm2d, status, ierr)
      call MPI_SENDRECV(cv(sx,ey+1),nx,MPI_DOUBLE_PRECISION,nbrright,5,
     1                 cv(sx,sy),nx,MPI_DOUBLE_PRECISION,nbrleft,5,
     2                 comm2d, status, ierr)
      call MPI_SENDRECV(h(sx,sy),nx,MPI_DOUBLE_PRECISION,nbrleft, 6,
     1                 h(sx,ey+1),nx,MPI_DOUBLE_PRECISION,nbrright,6,
     2                 comm2d, status, ierr)
      call MPI_SENDRECV(z(sx+1,ey+1),nx,MPI_DOUBLE_PRECISION,nbrright,7,
     1                 z(sx+1,sy),nx,MPI_DOUBLE_PRECISION,nbrleft,7,
     2                 comm2d, status, ierr)     
*     Special comms between corner processes
*     Top right gets from bottom left
*     Bottom left gets from top right
*     Top left gets from bottom right
*     Bottom right gets from top left

* Top Right     
      if (myid .eq. topRight .and. myid .ne. bottomLeft) then
        call MPI_RECV(cu(sx,ey+1),1, MPI_DOUBLE_PRECISION, bottomLeft,
     1                 m,comm2d,status,ierr)
        call MPI_SEND(cv(sx,ey+1),1, MPI_DOUBLE_PRECISION,bottomLeft,
     1                 m+1,comm2d,ierr)
      end if
* Bottom Left     
      if (myid .eq. bottomLeft.and. myid .ne. topRight) then
        call MPI_SEND(cu(ex+1,sy),1, MPI_DOUBLE_PRECISION,topRight,
     1                 m,comm2d,ierr)
        call MPI_RECV(cv(ex+1,sy),1, MPI_DOUBLE_PRECISION, topRight,
     1                 m+1,comm2d,status,ierr)
      end if
* Single processor case     
      if (myid .eq. bottomLeft .and. myid .eq. topRight) then
       CU(1,N+1) = CU(M+1,1)
       CV(M+1,1) = CV(1,N+1)
      end if     
*Top Left     
      if (myid .eq. topLeft .and. myid .ne. bottomRight) then
         call MPI_RECV(z(sx,sy),1, MPI_DOUBLE_PRECISION, bottomRight,
     1                 m,comm2d,status,ierr)
         call MPI_SEND(h(sx,sy),1, MPI_DOUBLE_PRECISION,bottomRight,
     1                 m+1,comm2d,ierr)
      end if
*Bottom Right     
      if (myid .eq. bottomRight .and. myid .ne. topLeft) then
         call MPI_SEND(z(ex+1,ey+1),1, MPI_DOUBLE_PRECISION,topLeft,
     1                 m,comm2d,ierr)
         call MPI_RECV(h(ex+1,ey+1),1, MPI_DOUBLE_PRECISION, topLeft,
     1                 m+1,comm2d,status,ierr)
      end if
      if (myid .eq. bottomRight .and. myid .eq. topLeft) then
       Z(1,1) = Z(M+1,N+1)
       H(M+1,N+1) = H(1,1)
      end if
      end subroutine capuvhz

      END
C


c
c This routine shows how to determine the neighbors in a 2-d decomposition of
c the domain. This assumes that MPI_Cart_create has already been called 
c
      subroutine fnd2dnbrs( comm2d, 
     $                      nbrleft, nbrright, nbrtop, nbrbottom )
      integer comm2d, nbrleft, nbrright, nbrtop, nbrbottom
c
      integer ierr
c
      call MPI_Cart_shift( comm2d, 1,  1, nbrleft,   nbrright, ierr )
      call MPI_Cart_shift( comm2d, 0,  1, nbrtop, nbrbottom,   ierr )
c  This was in original code - did I miss something in translation??
*     call MPI_Cart_shift( comm2d, 0,  1, nbrleft,   nbrright, ierr )
*     call MPI_Cart_shift( comm2d, 1,  1, nbrbottom, nbrtop,   ierr )
c
      return
      end
c
      subroutine fnd2ddecomp( comm2d, m, n, sx, ex, sy, ey )
      integer comm2d
      integer m, n, sx, ex, sy, ey
      integer dims(2), coords(2), ierr
      logical periods(2)

c
      call MPI_Cart_get( comm2d, 2, dims, periods, coords, ierr )

      call MPE_DECOMP1D( m, dims(1), coords(1), sx, ex )
      call MPE_DECOMP1D( n, dims(2), coords(2), sy, ey )
c
      return
      end
c
c  This file contains a routine for producing a decomposition of a 1-d array
c  when given a number of processors.  It may be used in "direct" product
c  decomposition.  The values returned assume a "global" domain in [1:n]
c
      subroutine MPE_DECOMP1D( n, numprocs, myid, s, e )
      integer n, numprocs, myid, s, e
      integer nlocal
      integer deficit
c
      nlocal  = n / numprocs
      s       = myid * nlocal + 1
      deficit = mod(n,numprocs)
      s       = s + min(myid,deficit)
      if (myid .lt. deficit) then
          nlocal = nlocal + 1
      endif
      e = s + nlocal - 1
      if (e .gt. n .or. myid .eq. numprocs-1) e = n
      return
      end

      character*4 function int_to_string(int)
*
* This needs the usual checks: 0 <= int < 10000
*     
      integer int
      integer digit, temp, i
      character(len=4):: out_string='0000'
      integer base
      base = iachar('0') 
      temp = int
      i = 4
      do while (temp > 0 .and. i > 0)
         digit = mod(temp,10)
         out_string(i:i) = achar(base+digit)
         temp = temp / 10
         i = i - 1
      end do
      int_to_string = out_string     
      end

