/*integral_trapezium.c*/

/*Integration using Trapezium rule*/


#include <stdio.h>
#include <math.h>
#include <float.h>

#define N 200
#define PI 3.1415927



double testfunction (double x)
{
    double y=sin(x);
    return y;
}




int main(int argc, char **argv)
{
	int i,n;
	double a,b,h,sum,result;
	double y[N];
        double x[N];


	/*integration interval*/
	a=0;
	b=PI/2;

	h=(b-a)/(N+1);


	for(i=0;i<N;i++)
	{
	    x[i]=a+h*i;
	    y[i]=testfunction(x[i]);
	}

	sum=0;


	//quadrature using the trapezoid rule
	//accuracy to order h^2
	for(i=0;i<=N-2;i++)
	    sum=sum+y[i]+y[i+1];
	

	result=sum*h/2;

	printf("results is %g\n",result);

}
