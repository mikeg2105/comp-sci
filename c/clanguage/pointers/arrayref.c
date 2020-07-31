/*program12.c Different methods for referencing array elements*/
/*arrayref.c*/
#include <stdio.h>

main()
{
	int i, offset;
        int iVec[] = {11,22,33,44,55};
	int *iVecptr = iVec;  /*Set iVecptr to the iVec array*/

	printf("Displaying Array iVec using array subscripts:\n");
        for(i=0; i<4; i++)
          printf("iVec[%d]=%d\n", i, iVec[i]);

	printf("Displaying Array iVec using pointer offsets:\n");
        for(offset=0; offset<4; offset++)
          /*Note: Use of * operator to dereference the
            offset pointer to the array*/
          printf("*(iVec+%d)=%d\n", offset, *(iVec+offset));

	printf("Displaying Array iVec using pointer subscripts:\n");
        for(i=0; i<4; i++)
          printf("iVecptr[%d]=%d\n", i, iVecptr[i]);

	printf("Displaying Array iVec using pointer offsets:\n");
        for(offset=0; offset<4; offset++)
          /*Note: Use of * operator to dereference the
            offset pointer to the array*/
          printf("*(iVecptr+%d)=%d\n", offset, *(iVecptr+offset));

	return 0;
}
