/*program9.c Examples of function calls */
/*           Call by reference          */
/*           Call by value              */	

/*functioncalls.c*/
#include <stdio.h>

float FuncByValue(float);
void  FuncByReference(float *);

main()
{
	float finval = 4;
        float finref = 5;

        printf("Original value of finval=%f\n", finval);
        finval=FuncByValue(finval);
        printf("The new value of finval after by value call is %f\n", finval);

        printf("Original value of finref=%f\n", finref);
        FuncByReference(&finref);
        printf("The new value of finref after by value call is %f\n", finref);

	return 0;
}

float FuncByValue(float fval)
{
	return fval*fval;
}

void  FuncByReference(float *fvalptr)
{
	*fvalptr = *fvalptr * *fvalptr;
}
