/*integral_simpson.c*/

/*Integration using Simpsons rule*/


#include <stdio.h>
#include <math.h>
#include <float.h>

#define N 200
#define PI 3.1415927

//quadrature to higher order using lagrange interpolation

double testfunction (double x)
{
    double y=sin(x);
    return y;
}

//lagrange interpolation

//xval is the value of x at which we wish to interpolate
//f is the series of values
double lagrange_interp(double xval,double *f, double *x, int i)
{  
    double y,t1,t2,t3;
    t1=(xval-x[i])*(xval-x[i+1])/((x[i-1]-x[i])*(x[i-1]-x[i+1]));
    t2=(xval-x[i-1])*(xval-x[i+1])/((x[i]-x[i-1])*(x[i]-x[i+1]));
    t3=(xval-x[i-1])*(xval-x[i])/((x[i+1]-x[i-1])*(x[i+1]-x[i]));
    y=t1*f[i-1]+t2*f[i]+t3*f[i+1];
    return y;
}


int main(int argc, char **argv)
{
	int i,n;
	double a,b,h,sum,result;
	double y[2*N+1];
        double x[2*N+1];

	/*use an odd number of elements*/
	n=(N*2)+1;

	/*integration interval*/
	a=0;
	b=PI/2;

	h=(b-a)/n;


	for(i=0;i<n;i++)
	{
	    x[i]=a+h*i;
	    y[i]=testfunction(x[i]);
	}

	sum=0;

	//since we are using interpolation
	//the loop performs the summation 
	//over every other element
	for(i=1;i<=n-1;i+=2)
	    sum=sum+y[i-1]+4*y[i]+y[i+1];
	

	result=sum*h/3;

	printf("results is %g\n",result);

}
