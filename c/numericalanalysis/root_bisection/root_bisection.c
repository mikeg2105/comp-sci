/*Bisection method for finding roots*/

/* here is an example use of the while statement
// which is used for finding the root of a polynomial 
// which is known to lie within a certain interval.
// a is the lower value of the range
// b is the upper value of the range */
#include <stdio.h>
#include <math.h>
#include <float.h>


/*
  Note FLT_MIN, FLT_MAX and FLT_EPSILON
  defined in float.h
*/


float sign(float f){return(fabs(f)/f);}

int main(char **argv, int argc)
{

	float x,fx;
	float a = 0;
	float fa = -FLT_MIN;
	float b = 3;
	float fb = FLT_MAX;

	while( fabs(b-a)>(FLT_EPSILON*b))
	{

	  x = (a+b)/2;

	  /*The function whose root is to be determined*/
	  fx = pow(x,3)-2*x-5;
	  if(sign(fx)==sign(fa))
	  {
	    a = x;
	    fa = fx;
            printf("a=%f fa=%f\n",a,fa);
	  }
	  else
	  {
	    b = x;
	    fb = fx;
            printf("b=%f fb=%f\n",b,fb);
	  }

	}

	printf(" The root is :%f\n",x);

	return 0;
}
