/*program10: Initialisin and using an array*/
/*array.c*/


#include <stdio.h>

#define VECSIZE 5
#define MATROW 3
#define MATCOL 3


float SumVec(float []);
float SumMatDiagonal(float[][]);
 
main()
{
	int i, j; /*Loop control variables*/
        float fvecsum, fmattrace;     /*Output data*/
        	
	float myvec[VECSIZE] = {1,2,3,4,5};

	float mymat[MATROW] [MATCOL]={{1,2,4},
				      {3,5,8},
                                      {5,8,12}};

	printf("Your input vector is:\n");
        for(i=0; i<VECSIZE; i++)
           printf("%f ", myvec[i]);
        printf("\n");
        fvecsum=SumVec(myvec);
        printf("The vector sum is: %f\n\n", fvecsum);

	printf("Your input matrix is:\n");
        for(i=0; i<MATROW; i++)
	{
		for(j=0; j<MATCOL; j++)
			    printf("%f ", mymat[i][j]);
		printf("\n");
	}
        printf("\n");
        fmattrace=SumMatDiagonal(mymat);
	printf("The trace of the matrix is %f\n", fmattrace);
    
	return 0;
}

float SumVec(float fVec [])
{
	float sum=0;
        int i;

	for(i=0; i<VECSIZE; i++)
		sum += fVec [i];

	return sum;
}

float SumMatDiagonal(float fMat[][MATCOL])  /*Note the number of columns is */
					    /*  defined here*/
{
	float diagonal=0;
        int i,j;

	for(i=0; i<MATROW; i++)
	{
                if(i>=MATCOL) break;  /*Break out of the for loop*/
                                          /*if the column index exceeds MATCOL*/
		diagonal = diagonal + fMat[i][i];
	}

	return diagonal;
}
