/*program11.c */
/*Use bubble sort to sort an
  array in descending order*/
  
 /*bubblesort.c*/

#include <stdio.h>

#define ARRAYSIZE 15


void BubbleSort(float *);
void Swap(float *, float *);
 
main()
{
	int i; /*Loop control variables*/
        	
	float myarray[ARRAYSIZE] = {3,23,552,6395,2828,54,5,12,32,64,11,55,3425,143,673};


	printf("Before sorting my array is:\n");
        for(i=0; i<ARRAYSIZE; i++)
           printf("%f ", myarray[i]);
        printf("\n\n");

        BubbleSort(myarray);

	printf("After sorting my array is:\n");
        for(i=0; i<ARRAYSIZE; i++)
           printf("%f ", myarray[i]);
        printf("\n");

	return 0;
}

void BubbleSort(float *fArray)
{
	int pass, j;

	for(pass=1; pass<=ARRAYSIZE-1; pass++)
			for(j=0; j<=ARRAYSIZE-2;j++)
				if(fArray[j]<fArray[j+1])
					Swap(&fArray[j], &fArray[j+1]);

}

void Swap(float *f1ptr, float *f2ptr)  
{
	float temp;

	temp = *f1ptr;
	*f1ptr = *f2ptr;
        *f2ptr = temp;
}

