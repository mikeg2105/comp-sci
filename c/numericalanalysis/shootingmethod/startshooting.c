
/*PROGRAM STARTSHOOTING
  THIS PROVIDES A SOLUTION FOR A BOUNDARY-VALUE PROBLEM FOR A
  FIRST-ORDER ODE WITH ONE UNKNOWN PARAMETER TO BE DETERMINED.
  THE ODE IS OF THE FORM     Dy/Dx = f(x,Q)  WHERE Q IS THE 
  UNKNOWN PARAMETER.   THE BOUNDARY CONDITIONS ARE y = Ya FOR
  x = a AND y = Yb FOR X = b.   
  THE FUNCTION f(x,Q) IS PROVIDED AS A FUNCTION STATEMENT.*/
#include <stdio.h>
#include <math.h>


double fun(double x,double q)
{
	return (-15.915494*q/pow((2-x),2));
}

int main(int argc, char **argv)
{
	double y[1000],dy[4];
        double x,xx,xi,xf,yf,q,h;
        double w[4]={0.0,0.5,0.5,1.0};
	char ans;
	FILE *fd;

	char instring[256];

        int i,j,n,k,m;

        double *qptr=&q;

        xi=0;
        y[0]=500;

        printf("THE FIRST BOUNDARY CONDITION  \n");
        printf(" y = %g WHEN x = %g.\n",y[0],xi);

	xf=1.0;
	yf=300;

        printf("THE SECOND BOUNDARY CONDITION \n");
        printf("WHERE y = %g WHEN x = %g.\n",yf,xf);


	n=100;
	m=20;
        h=(xf-xi)/n;


        //ESTIMATE OF THE UNKNOWN PARAMETER, Q.
        q=30.0;  //too high
        //q=20;   //too low     
        //q=25.1325;

        printf("Input your estimate for the initial value\n");
        scanf("%lf",&q);
        printf("Your q value is %lf \n",q);



        /*THE RUNGE-KUTTA INTEGRATION NOW BEGINS.F*/
        for(i=1;i<n;i++)
        {
		x=(i-1)*h;
		for(j=0;j<4;j++)
		{
			xx=x+h*w[j];
			
                        dy[j]=h*fun(xx,q);
		}
		y[i]=y[i-1]+(dy[0]+dy[3]+2*(dy[1]+dy[2]))/6.0;
	    }

        m=5; //output interval
        k=n/m;
        fd=fopen("mydoc.txt","w");
        for(i=0; i<n;i+=k)
        {
		x=i*h;
                fprintf(fd,"%f %f\n",x,y[i]);
                printf("%f %f\n",x,y[i]);
        }
        fclose(fd);

	return 0;
}

