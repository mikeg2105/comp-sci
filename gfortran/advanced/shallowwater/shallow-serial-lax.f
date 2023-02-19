      	PROGRAM SHALLOW
      	use omp_lib
      	implicit none
C ***
C *** 	BENCHMARK WEATHER PREDICTION PROGRAM FOR COMPARING THE PERFORMANCE
C *** 	OF CURRENT SUPERCOMPUTERS. THE MODEL IS BASED ON THE PAPER, "THE
C *** 	DYNAMICS OF FINITE-DIFFERENCE MODELS OF THE SHALLOW-WATER
C *** 	EQUATIONS," BY R. SADOURNY, J. ATM. SCI., VOL.32, NO.4, APRIL 1975
C *** 	CODE BY PAUL N. SWARZTRAUBER, NATIONAL CENTER FOR ATMOSPHERIC
C *** 	RESEARCH, BOULDER, COLORADO   OCTOBER 1984.
*
C 	Steven McHale
C 	Tsunami Model
C 	Shallow-Water Wave Equation
C 	Crank-Nicholson Discretization

      	double precision g, u0, v0, h0, wavespeed
      	integer ni, nj, nt
      	double precision xmax, dx, ymax, dy
      	double precision dt, tmax, courant
      	double precision,allocatable :: u(:,:,:),v(:,:,:),b(:,:),h(:,:,:)
      	double precision,allocatable :: t(:),x(:),y(:)


C 	Constants
        g  = 9.81;
        u0 = 0;
        v0 = 0;
        b  = 0;
        h0 = 5030;

C 	Define the x domain
        ni = 151;
        xmax = 100000;
        dx = xmax/(ni-1);

C 	Define the y domain
        nj = 151;
        ymax = 100000;
        dy = ymax/(nj-1);

C 	Define the wavespeed
        wavespeed = u0 + sqrt(g*(h0 - b));

C 	Define time-domain
        dt = 0.68*dx/wavespeed;
        tmax = 1500;

        nt=tmax/dt
        courant = wavespeed*dt/dx;

      	allocate (u(1:ni,1:nj,1:nt))
      	allocate (v(1:ni,1:nj,1:nt))
      	allocate (h(1:ni,1:nj,1:nt))
      	allocate (b(1:ni,1:nj))
      	allocate (t(1:nt))
      	allocate (x(1:ni))
      	allocate (y(1:nj))

	do i = 1 , ni
!		x(1:ni) = dx*[i,i=1,i=ni]
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


C Build empty u, v, b matrices
!    u(1:ni,1:nj,1:nt)=0
!    v(1:ni,1:nj,1:nt)=0
!    h(1:ni,1:nj,1:nt)=0
!    b(1:ni,1:nj)=0







    endprogram



