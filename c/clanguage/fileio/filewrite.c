/*program20.c*/
/*Example of Writing data to a sequential file */

/*filewrite.c*/

#include <stdio.h>
#include <stdlib.h>


main()
{
	int iRunid;
	char sName[30];
        int iNode;

	FILE *mfptr; /* mfptr = models.dat file pointer*/

	printf("Simulation log generator.\n");

       
        if((mfptr=fopen("models.dat", "w"))==NULL)
		printf("The file could not be opened!\n");
	else
	{
		printf("Enter the runid, simulation name, and node id.\n");
		printf("Enter EOF toterminate input.\n");
		printf("? ");
		scanf("%d%s%d", &iRunid, sName, &iNode);
		while(!feof(stdin))
		{
			fprintf(mfptr, "%6d %20s %6d\n", iRunid, sName, iNode);
			printf("? ");
			scanf("%d%s%d", &iRunid, sName, &iNode);
		}

		fclose(mfptr);
	}

	return 0;
	

}
