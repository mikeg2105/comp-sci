//Determne the root of a function using the Newton-Raphson method
// used in association with the finding of a root
// problem solved earlier.
// % a is the lower value of the range
// b is the upper value of the range

#include <stdio.h>
#include <math.h>
#include <float.h>

#define PI 3.1415927


int main(char **argv, int argc)
{
	float x=3;
        float fx;
	float a = 0;
	float fa = -FLT_MIN;
	float b = 3;
	float fb = FLT_MAX;
	float dfx;

	printf("Computing a root using the Newton-Raphson method\n");
	while(   (  sqrt((x-a)*(x-a))   )>FLT_EPSILON)
	{
	  a=x;
	  //fx = x^3-2*x-5;
	  //dfx= 3*x^2-2;
	  
	  fx = 0.5*sin(2*(x-(PI/4)))+0.5*sin(x);
	  dfx= cos(2*(x-(PI/4)))+0.5*cos(x);
	  x=a-(fx/dfx);
	  printf("x-a=%g\n",x-a);
	} 

	printf(" The root is : %g\n",x);

	return 0;
}
