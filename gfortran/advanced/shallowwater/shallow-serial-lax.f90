      	PROGRAM SHALLOW
      	use omp_lib
      	implicit none
!! ***
! *** 	BENCHMARK WEATHER PREDICTION PROGRAM FOR COMPARING THE PERFORMANCE
! *** 	OF CURRENT SUPERCOMPUTERS. THE MODEL IS BASED ON THE PAPER, "THE
! *** 	DYNAMICS OF FINITE-DIFFERENCE MODELS OF THE SHALLOW-WATER
! *** 	EQUATIONS," BY R. SADOURNY, J. ATM. SCI., VOL.32, NO.4, APRIL 1975
! *** 	CODE BY PAUL N. SWARZTRAUBER, NATIONAL CENTER FOR ATMOSPHERIC
! *** 	RESEARCH, BOULDER, COLORADO   OCTOBER 1984.
!
! 	Steven McHale
! 	Tsunami Model
! 	Shallow-Water Wave Equation
! 	Crank-Nicholson Discretization

      	double precision g, u0, v0, h0, b0, wavespeed
      	doubleprecision loweri, lowerj
      	integer ni, nj, nt
      	integer i,j,n
      	double precision xmax, dx, ymax, dy
      	double precision dt, tmax, courant
      	double precision,allocatable :: u(:,:,:),v(:,:,:),b(:,:),h(:,:,:)
      	double precision,allocatable :: t(:),x(:),y(:)

        character(len=4) int_to_string
! 	Constants
        g  = 9.81
        u0 = 0
        v0 = 0
        b0  = 0
        h0 = 5030;

! 	Define the x domain
        ni = 151
        xmax = 100000
        dx = xmax/(ni-1)

! 	Define the y domain
        nj = 151
        ymax = 100000
        dy = ymax/(nj-1)

! 	Define the wavespeed
        wavespeed = u0 + sqrt(g*(h0 - b0))

! 	Define time-domain
        dt = 0.68*dx/wavespeed
        tmax = 1500

        nt=tmax/dt
        courant = wavespeed*dt/dx

      	allocate (u(1:ni,1:nj,1:nt))
      	allocate (v(1:ni,1:nj,1:nt))
      	allocate (h(1:ni,1:nj,1:nt))
      	allocate (b(1:ni,1:nj))
      	allocate (t(1:nt))
      	allocate (x(1:ni))
      	allocate (y(1:nj))

	do i = 1 , ni
!		x(1:ni) = dx*(i,i=1,i=ni)
		x(i) = dx*i
	enddo

	do j=1,nj
!		y(1:nj) = dy*[j,j=1,j=nj]
		y(j) = dy*j
	enddo

	do i=1,nt
!		t(1:nt) = dt*[i,i=1,i=nt]
		t(i) = dt*i
	enddo


! Build empty u, v, b matrices
    u(1:ni,1:nj,1:nt)=0
    v(1:ni,1:nj,1:nt)=0
    h(1:ni,1:nj,1:nt)=0
    b(1:ni,1:nj)=0

!Define h
    loweri=55000/100000*(ni-1)+1
    lowerj=55000/100000*(nj-1)+1
    h(:,:,1) = 5000;
!    h((45000/100000*(ni-1)+1):floor(loweri),(45000/100000*(nj-1)+1):floor(lowerj),1) = 5030;
    h(68:83,68:83,1) = 5030;
!Define b
    do i = 1,ni
        if (x(i) > 20001) then
            b(:,i) = 0;
        else if (x(i) < 20000) then
            b(:,i) = 5000/20000*(20000-x(i));
        end if
    enddo

    print *, "initilisation complete"
    print *,(45000/100000*(ni-1)+1),floor(loweri),(45000/100000*(nj-1)+1),floor(lowerj)
!     Employ Lax
    do n=1,nt-1
        do i=2,(ni-1)
            do j=2,(nj-1)
                u(i,j,n+1) = ( ( u(i+1,j,n) + u(i-1,j,n) + u(i,j+1,n) + u(i,j-1,n) )/4 ) &
                 - 0.5*(dt/dx)*( ( u(i+1,j,n)**2)/2 - (u(i-1,j,n)**2)/2 ) &
                 - 0.5*(dt/dy)*( v(i,j,n)) *( u(i,j+1,n) - u(i,j-1,n)) - 0.5*g*(dt/dx)*( h(i+1,j,n)-h(i-1,j,n))

                v(i,j,n+1) = ((v(i+1,j,n) + v(i-1,j,n) + v(i,j+1,n) + v(i,j-1,n))/4) &
                - 0.5*(dt/dy)*((v(i,j+1,n)**2)/2 - (v(i,j+1,n)**2)/2) &
                - 0.5*(dt/dx)*(u(i,j,n))*(v(i+1,j,n) - v(i-1,j,n)) - 0.5*g*(dt/dy)*(h(i,j+1,n)-h(i,j-1,n))

                h(i,j,n+1) = ((h(i+1,j,n) + h(i-1,j,n) + h(i,j+1,n) + h(i,j-1,n))/4) &
                - 0.5*(dt/dx)*(u(i,j,n))*((h(i+1,j,n)-b(i+1,j)) - (h(i-1,j,n)-b(i-1,j))) &
                - 0.5*(dt/dy)*(v(i,j,n))*((h(i,j+1,n)-b(i,j+1)) - (h(i,j-1,n)-b(i,j-1))) &
                - 0.5*(dt/dx)*(h(i,j,n)-b(i,j))*(u(i+1,j,n)- u(i-1,j,n)) &
                - 0.5*(dt/dy)*(h(i,j,n)-b(i,j))*(v(i,j+1,n) - v(i,j-1,n))

            enddo
        enddo

!     Define Boundary Conditions
    u(1,:,n+1) = 2.5*u(2,:,n+1) - 2*u(3,:,n+1) + 0.5*u(4,:,n+1)
    u(ni,:,n+1) = 2.5*u(ni-1,:,n+1) - 2*u(ni-2,:,n+1) + 0.5*u(ni-3,:,n+1)
    u(:,1,n+1) = 2.5*u(:,2,n+1) - 2*u(:,3,n+1) + 0.5*u(:,4,n+1)
    u(:,nj,n+1) = 2.5*u(:,nj-1,n+1) - 2*u(:,nj-2,n+1) + 0.5*u(:,nj-3,n+1)

    v(1,:,n+1) = 2.5*v(2,:,n+1) - 2*v(3,:,n+1) + 0.5*v(4,:,n+1);
    v(ni,:,n+1) = 2.5*v(ni-1,:,n+1) - 2*v(ni-2,:,n+1) + 0.5*v(ni-3,:,n+1);
    v(:,1,n+1) = 2.5*v(:,2,n+1) - 2*v(:,3,n+1) + 0.5*v(:,4,n+1);
    v(:,nj,n+1) = 2.5*v(:,nj-1,n+1) - 2*v(:,nj-2,n+1) + 0.5*v(:,nj-3,n+1);

    h(1,:,n+1) = 2.5*h(2,:,n+1) - 2*h(3,:,n+1) + 0.5*h(4,:,n+1);
    h(ni,:,n+1) = 2.5*h(ni-1,:,n+1) - 2*h(ni-2,:,n+1) + 0.5*h(ni-3,:,n+1);
    h(:,1,n+1) = 2.5*h(:,2,n+1) - 2*h(:,3,n+1) + 0.5*h(:,4,n+1);
    h(:,nj,n+1) = 2.5*h(:,nj-1,n+1) - 2*h(:,nj-2,n+1) + 0.5*h(:,nj-3,n+1);


         open(1,form='unformatted', &
               file=trim(".")// '/'//trim("data//test")//int_to_string(n)//trim(".out"))
         write(1) h(1:ni,1:nj,n)
         close(1)
    enddo




    endprogram

    character*4 function int_to_string(int)
!
! This needs the usual checks: 0 <= int < 10000
!
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

