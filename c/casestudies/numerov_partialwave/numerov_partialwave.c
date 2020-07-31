
//An illustration fo the numerov method for solving partial diffrential equations
//In this case the Schrodinger equation is solved.
//The problem here is the use of partial wave analysis to model the scattering of
//a particle interacting with a potential well.

//In this case the potential well is the Lennard-Jones potential used to represent 
//atomic interactions for example between inert atoms such as Argon



//http://en.wikipedia.org/wiki/Numerov's_method

#include <stdio.h>
#include <math.h>
#include <float.h>

#define N 100
#define PI 3.1415927

#ifdef USE_GSL
        /*
          gcc -Wall -I/usr/local/include -c example.c
          gcc -L/usr/local/lib example.o -lgsl -lgslcblas -lm
	*/
	#include <gsl/gsl_math.h>
	#include <gsl/gsl_sf_bessel.h>
#endif
#ifdef USE_NAG
	#include <nag.h>
	#include <nag_stdlib.h>
	#include <nags.h>
#endif


//Global variables used in this program 
//and by all functions.  Naughty.....?

double m;   //mass of particle being scattered
double hb; //Plancks constant over 2*pi


/*Based on taylor expansion about 0 i.e. valid for mod(x)<1 */

#ifdef USE_NAG
double nagbessel(int l, double x)
{
	double bes;
switch(l)
{
	case 0:
		bes=s17aec(x, NAGERR_DEFAULT);
	break;
	case 1:
		bes=s17afc(x, NAGERR_DEFAULT);
	break;
	case 2:
		bes=(pow(x,2.0)/8)+(pow(x,4.0)/96)+(pow(x,6.0)/3072);
	break;
	case 3:
		bes=(pow(x,3.0)/48)+(pow(x,5.0)/768);
	break;
	case 4:
		bes=(pow(x,4.0)/384)+(pow(x,6.0)/7680);
	break;

	default:
		bes=(pow(x,5.0)/3840);
	break;
}

	return bes;
}
#endif

/*Based on taylor expansion about 0 i.e. valid for mod(x)<1 */
double besselapprox(int l, double x)
{
	double bes;
switch(l)
{
	case 0:
		bes=1+(pow(x,2.0)/4)+(pow(x,4.0)/64)+(pow(x,6.0)/2304);
	break;
	case 1:
		bes=(x/2)+(pow(x,3.0)/16)+(pow(x,5.0)/384);
	break;
	case 2:
		bes=(pow(x,2.0)/8)+(pow(x,4.0)/96)+(pow(x,6.0)/3072);
	break;
	case 3:
		bes=(pow(x,3.0)/48)+(pow(x,5.0)/768);
	break;
	case 4:
		bes=(pow(x,4.0)/384)+(pow(x,6.0)/7680);
	break;

	default:
		bes=(pow(x,5.0)/3840);
	break;
}

	return bes;
}




double tdl(double u1,double u2,double r1,double r2,int l,double k)
{
 //calculate phase shift of partial waves tand deltl
  double kk=(r1*u1)/(r2*u2);

#ifdef USE_GSL
 double tdlv=(kk*gsl_sf_bessel_Jn(l,k*r1)-gsl_sf_bessel_Jn(l,k*r2))/(kk*gsl_sf_bessel_Jn(l,k*r1)-gsl_sf_bessel_Jn(l,k*r2));
#else
	#ifdef USE_NAG
	 double  tdlv=(kk*nagbessel(l,k*r1)-nagbessel(l,k*r2))/(kk*nagbessel(l,k*r1)-nagbessel(l,k*r2));

	#else
	 double  tdlv=(kk*besselapprox(l,k*r1)-besselapprox(l,k*r2))/(kk*besselapprox(l,k*r1)-besselapprox(l,k*r2));
	#endif
#endif
	return tdlv;
}




double u(int l, double r)
{
  double alpha=6.12;
  double e=3;
  double C=sqrt(e*alpha/25);
  double uv=exp( (-1)*C*pow(r,-5.0));
   return uv;
}




double v(double r, double sigma, double epsilon)
{
   //epsilon=5.9;//meV H-Kr interaction
  //sigma=3.57;//Angstrom
  double vv=10*epsilon*(    (  pow((sigma/r),12.0)  )-2*(  pow((sigma/r),6.0)  )       );
    return vv;
}





double f(int l, double r, double e,double sigma,double epsilon)
{
	double fv=v(r,sigma,epsilon)+((hb*hb)/(2*m*(r*r)))-e;
	return fv;
}


//solve radial wave equation using numerov 
//predictor method to fourth order 
//for 2nd order de's

double numerov(double u1,double u2,double l,double r,double delta,double e,double sigma,double epsilon)
{
//solve radial wave equation using numerov 
//predictor method to fourth order 
//for 2nd order de's

  double num1=1/(1-(pow(delta,2.0)/12)*f(l,r+delta,e,sigma,epsilon));
  double bracket1=2*u2-u1+((delta*delta)/12)*(10*f(l,r,e,sigma,epsilon)*u2+f(l,r,e,sigma,epsilon)*u1);
  double numerovv=num1*bracket1;

   return numerovv;
}




int main(int argc, char **argv)
{
double deltah=0.01;
int nsteps=200;
int i,j,nec;

FILE *mfptr; /* mfptr = file pointer for results file*/

//m=938*10^9;
m=1.672*pow(10.0,-27.0);
//hb=6.59*10^(-13);
hb=1.054*pow(10.0,-34.0);
int lupper=10;
double e=1.6*pow(10.0,-19.0);
 
double delta=0.1*pow(10.0,-10.0);
 int nr,ne;
double k,sumdelta,epsilon,sigma,cosecdelta2,res;
 delta=0.5;
 m=1;
 hb=1;
 e=1;

//2m/hb^2=6.12meV^-1(sigma)^-2
  lupper=10;
  sumdelta=0;
  nr=100;
  ne=250;


//lupper=4;
//nr=10;
//ne=15;
  
  double *sumouter=(double *)calloc(ne,sizeof(double));
  double totsum=0;
  //sigma=0;
  double **u1=(double **)calloc(ne,sizeof(double *));
  double **u2=(double **)calloc(ne,sizeof(double *));
  double **u3=(double **)calloc(ne,sizeof(double *));

       if((mfptr=fopen("partialwave.dat", "w"))==NULL)
       {
		printf("The file could not be opened!\n");
                return 1;
       }


			fprintf(mfptr, "Energy \t Cross-Section\n");

  for(i=0; i<ne; i++)
  {
    u1[i]=(double *)calloc(lupper+1,sizeof(double));
    u2[i]=(double *)calloc(lupper+1,sizeof(double));
    u3[i]=(double *)calloc(lupper+1,sizeof(double));
  }





  //outer loop integration over r
  epsilon=5.9;//meV H-Kr interaction
  sigma=3.57;//Angstrom

  for( nec=0; nec<ne; nec++)
  {
    e=(nec+1)*0.0005;
    k=sqrt(2*m*e)/hb;
    totsum=0;
  //sigma=0;
  
  for(i=0; i<ne; i++)
  {
	  sumouter[i]=0.0;
	  for(j=0; j<lupper+1; j++)
	  {
		u1[i][j]=0.0;
		u2[i][j]=0.0;
		u3[i][j]=0.0;
	  }
  }

  for(j=0;j<=nr;j++)
  {
    
    //u1=0;
    //inner loop summation over l
    sumdelta=0;
    for(i=0; i<=lupper; i++)
    {   

      if(j == 0)
      {

        u1[j][i]=.1;
        u2[j][i]=pow(delta,i);
      }
      else
      {
        u2[j][i]=u3[j-1][i];
        u1[j][i]=u2[j-1][i];
      }
      
      u3[j][i]=numerov(u1[j][i],u2[j][i],i,3.1+j*delta,delta,e,sigma,epsilon);
     } //for(i=0; i<=lupper; i++)
    
  }  //for j=1:nr,
  
 
     for( i=0;i<=lupper; i++)
     {   
       res=tdl(u2[nr-1][i],u2[nr][i],(nr-1)*delta,(nr)*delta,i,k);
       cosecdelta2=((1/(pow(res,2.0)))+1);
       sumdelta=sumdelta+(2*i+1)*(1/cosecdelta2); 
     }
    //end
    sumouter[nec]=((4*PI)/(pow(k,2.0)))*sumdelta;
    totsum=totsum+sumouter[nec];
 
    fprintf(mfptr, "%g %g\n", (nec+1)*0.0005, sumouter[nec]);
  }//for( int nec=1; nec<=ne; nec++)  
   
 

		
	

  fclose(mfptr);
  
  free(sumouter);
  free(u1);
  free(u2);
  free(u3);


}
