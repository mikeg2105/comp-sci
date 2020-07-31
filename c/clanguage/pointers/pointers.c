/*Program8: Example of the use of pointer operators */
/*pointers.c*/

#include <stdio.h>

main()
{
	float f1;
        float *f1ptr;  /* Declare a pointer variable to an integer*/

	f1=3.141;
	f1ptr=&f1; /*f1ptr is set to the address of f1*/

	printf("The address of f1 is %p\n"
               "The value of f1ptr is %p\n\n", &f1, f1ptr);

	printf("The value of f1 is %f\n"
	       "The pointer f1ptr points to the value %f\n\n", f1, *f1ptr);

	printf("Proving that * and & are complementary operators.\n"
	       "&*f1ptr = %p\n *&f1ptr = %p\n", &*f1ptr, *&f1ptr);

	return 0;
}
