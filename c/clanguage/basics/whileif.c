/*Example 3: Example of while repetition and  
  conditional statements*/
/*whileif.c*/
#include <stdio.h>

main()
{
    /*initializing variables in decalarations*/
    int files=1, ntitania_files=0, nmaxima_files=0,
        nother_files=0, result;

    /*Process 5 files with a counter-controlled loop*/
    while(files<=5)
    {
	printf("Enter file location (1=Titania, 2=Maxima): ");
        scanf("%d", &result);

	if(result==1)
           ntitania_files = ntitania_files+1;
        else if(result==2)
           nmaxima_files = nmaxima_files+1;
        else
           nother_files=nother_files+1;

        files=files+1;

     }/*End of while file processing loop*/

     printf("There are %d files on Titania\n", ntitania_files);
     printf("There are %d files on Maxima\n", nmaxima_files);

     if(nother_files>2)
	printf("Who knows where the files are??\n");


     return 0;
}
