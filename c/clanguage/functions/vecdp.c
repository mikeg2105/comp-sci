/*program25.c*/
/*Use the time function to calculate the time
  required for vector dot product multiplication */

/*vecdp.c*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>


/*Wall clock function use time of day function to 
  return current time in seconds*/
void hpcwall(double *retval);
int doubvectest(int size, int nsteps);

int main(int argc, char* argv[])
{
	int nvsize, nmults;
	
	nvsize = 100; /*Size of vector*/
	nmults = 10000000; /*Number of multiplications*/

	if(argc>2)
	{
	   nvsize=atoi(argv[1]);
	   nmults=atoi(argv[2]);
	}
	printf("A vector of size %d \n" 
               "Will be %d times multiplied using a dot product.\n", nvsize, nmults);

	return doubvectest(nvsize, nmults);	
}

int doubvectest(int size, int nsteps)
{
   int status = 0;

	int i, j, k;

	time_t startt, finisht; 
	double totalt;

	double *veca, *vecb;
	int test;

	double result=0;

	FILE *ofile;

	/*Create vectors*/
	veca = (double *)calloc(size, sizeof(double));
	vecb = (double *)calloc(size, sizeof(double));

	printf("Created Vectors\n");

	if(ofile = fopen("doubvecresults.txt", "w"))
	{
		/*Fill vectors*/
		for(i=0; i<size; i++)
		{
			test = 1+rand()%1000;
			veca [i] = 1.0/test;
			test = 1+rand()%1000;
			vecb [i] = 1.0/test;
		}
		printf("Filled Vectors\n");

		/*Calculate thetime required
	          for multiplication for different vector sizes*/
		fprintf(ofile, "Vector dot product timings for %d cycles", nsteps);
		for(k=1; k<size; k++)
		{
			startt = time(NULL);
			for(i=0; i<nsteps; i++)
				for(j=0; j<k; j++)
					result += (veca [j] * vecb [j]);

			finisht = time(NULL);
			totalt = (double)difftime(finisht,startt);

			printf("size=%d time=%f\n", k, totalt);
			fprintf(ofile, "vec size %d result %f total time %f \n", k, result, totalt);
		}
		fclose(ofile);
	}
	free(veca);
	free(vecb);

   return status;
}

void hpcwall(double *retval)
{
	static long zsec=0;
	static long zusec=0;
	double esec;
	
	struct timeval tp;
	struct timezone tzp;
	
	gettimeofday(&tp, &tzp);
	
	if(zsec==0) zsec=tp.tv_sec;
	if(zusec==0) zusec=tp.tv_usec;
	
	/*Usual  and correct return number of seconds*/
	*retval=(tp.tv_sec - zsec)+(tp.tv_usec-zusec)*0.000001;
	
	/*Un conventional return elapsed time in micro seconds*/
	/**retval=1000000*(tp.tv_sec - zsec)+(tp.tv_usec-zusec);*/
	
}
