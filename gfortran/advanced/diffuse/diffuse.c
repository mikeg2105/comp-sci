/*	 								     *
 *									     *
 * EXAMPLE INTRODUCING : POINT TO POINT COMMUNICATIONS....		     *
 *									     *
 * This program will simulate the mechanics of diffusion as follows:	     *
 * Given that the number of processors  we can utilise is <n> 		     *
 *                               divide the domain into <n> regions. 	     *
 *  Impose an initial distribution of particles onto these <n> regions       *
 *			( for example a normal probability distribution).    *
 * Simulate diffusion by allowing <f> percent of particles from each domain  *
 * during a time-step to escape to the neighbouring domains on either side   *
 *		       of it. 1-D only so 2 neighbours. 		     *
 * NOTES: Particles that escape through each end never return.		     *
 * Simulate for a given number of time steps (m)			     *
 * List or plot the intermediate and final distributions.                    */
 
#include <stdio.h>
#include <mpi.h>
#include <math.h>	
	
#define ndims 1
   	
  int main(int argc, char *argv[])	
 {
        float normpdf( float x ,float mu , float sigma);
	int buffer[10000] ;
	int irank , ileft , iright , nsize , ierr , i , ai ;
	int  itstep=0 ;
	MPI_Status status;
	float dparams[3] , difpar , std , yright ,yleft ;
	float ynew , x ,ytot , ainterval;	 
	float periods[ndims] , dims[ndims]; 
	
	const float  dom_leng = 1.0;
 
	 MPI_Init( &argc , &argv );

	 MPI_Comm_size(MPI_COMM_WORLD, &nsize);

	 MPI_Comm_rank(MPI_COMM_WORLD, &irank ) ; 
/*	 MPI_Buffer_attach(&buffer,1000);   */
	
	if ( nsize <= 1 ) { 
	  printf ( " At least 2 processors are needed to run this job! ");
	  MPI_Finalize();
	  return (-1) ;
	    }
/*
 * initialization and user parameters input.
 */
	 
/*    master thread will read the user input and pass it onto others */		
	if ( irank == 0 ) {
	    printf( " total number of processors:%d", nsize );
	    printf( "\n  total length of the domain is taken to be 1 units");
            printf(" \n i.e. length is normalised. the initial distribution is assumed to ");
            printf(" \n center in the middle of the domain (i.e. at 0.5 ) " );
   
            std=-0.1;
            while (std <= 0.0 ) 
	      {
	      printf("\n Enter standard deviation of the initial distribution:" );
	      scanf("%f" ,  &std );
	      if( std <=0 ) printf("\nmust be positive"); 
	      } 

  
            printf("\n\nThe diffusion parameter will be specified as fraction of ");
	    printf("\nmass per unit length per time-step, must be between 0.0 and 1.0\n"); 
	    printf("\n Enter Diffusion Parameter:" ); 
	    scanf("%f" , &difpar ) ;

            difpar = difpar*dom_leng/nsize;		    
	    dparams[0] = std ;
	    dparams[1] = difpar	;            
        }
  
/* for good book-keeping: synchronise everything after reading data.  */  

	 MPI_Barrier( MPI_COMM_WORLD );
         

/* the input parameters ( dparams array) has only been read into thread 0 .  
   now pass them onto all the other threads via send and receives.   */
 	
	if ( irank == 0 ) 
 	  		for ( i = 1 ; i <= nsize-1 ; i++ ) 
                  MPI_Ssend( &dparams , 2 , MPI_FLOAT , i , 0 , MPI_COMM_WORLD); 	 
    else
    {
	      MPI_Recv(  &dparams , 2 , MPI_FLOAT , 0 , 0 ,MPI_COMM_WORLD,&status);
	      std = dparams[0];
	      difpar = dparams[1];
	}
 
	   
/* for good book-keeping: synchronise everything after reading data. */

	 MPI_Barrier( MPI_COMM_WORLD );
	  printf("\n std and difpar for rank %d = %f  %f",irank,std,difpar);

/* divide the domain amongst processes .	    
 ! we will use the mid-point of each interval in calculations.
 ! starting with ainterval/2 as the first intervals location. */
 
	 ainterval = dom_leng/ nsize ; 
	 		
/*	    	
		   
! distribute the mass amongst the processors.
!
!             --
!           --  -- 
!         --      ---
!      ---           ---
!  ----                 ------
!  -----------------------------------------      
! 1 , 2  , 3 , .....             n                       /*



/* initialise distribution during the first time step.  */

          x = ainterval*(irank + 0.5 ) ;
	  ytot = normpdf( x , dom_leng*0.5  , std ) ;
	  printf("\nFor rank %d the initial distribution at x= %f is %f",irank, x , ytot );   
	 
/* a time-steps.                  */ 

      for ( itstep = 1 ; itstep <= 10 ; itstep++ ) 
      {
		  /* by setting neighbours of end cells to null_proc we are assuming 
			 that particles migrate out as usual but never come back..            */
 
	      ileft = irank -1 ;
		  if (ileft < 0 ) ileft = MPI_PROC_NULL ;
		  iright = irank + 1 ;
		  if (iright > (nsize-1) ) iright = MPI_PROC_NULL ;
		  MPI_Sendrecv( &ytot ,1 ,MPI_FLOAT,iright,0,&yright,1,MPI_FLOAT,iright,0,
		                 MPI_COMM_WORLD,&status );
		   
		  MPI_Sendrecv( &ytot ,1 ,MPI_FLOAT,ileft,0,&yleft,1,MPI_FLOAT,ileft,0,
		                 MPI_COMM_WORLD, &status) ; 


	  ynew = ytot*(1.0-2.0*difpar) + difpar*(yleft+yright) ;
	  ytot = ynew ;
	  printf("\nTime Step %10d  Rank: %d   Distribution: %8f" ,itstep , irank , ytot );  

/*
 * the call below stops treads running away from each other 
 * commment it out and observe the effect ...                 */
 
	   MPI_Barrier( MPI_COMM_WORLD);	 
       }                               
      MPI_Finalize(); 
      printf("\n" ); 
      return (0) ;
     }
  
     float normpdf( float x ,float mu , float sigma) 
	 {
	   double  aaa, pi , asigma ;
	   pi= 4.0*atan(  (double) 1.0) ;
	   printf( "\n pi = %f " ,pi ) ;
	   asigma = (double) sigma;
	   if (asigma <= 0.0 ) asigma = 0.001;
	   aaa = -0.5* pow ( ((x-mu)/asigma) , (double)2.0 ) ;
	   printf(" partial %f" , exp(aaa) );
	   aaa = exp(aaa )/( asigma*sqrt(2.0*pi) ) ;
	   return ( (float) aaa) ;
	} 	
   
