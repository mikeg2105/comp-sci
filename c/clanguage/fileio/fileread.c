/*program21.c*/
/*Example of Reading data from a sequential file */
/*fileread.c*/

#include <stdio.h>

struct siminfo
{
	int id;
	char name[30];
	int node;
};

typedef struct siminfo SIMINFO;

main()
{

	SIMINFO sim;

	FILE *mfptr; /* mfptr = models.dat file pointer*/

	printf("Simulation log reader.\n");

       
        if((mfptr=fopen("models.dat", "r"))==NULL)
		printf("The file could not be opened!\n");
	else
	{
		printf("%6s %20s %6s\n", "Runid", "Name", "Node");
		fscanf(mfptr, "%d%s%d", &sim.id, sim.name, &sim.node);
		while(!feof(mfptr))
		{
			printf("%6d %20s %6d\n", sim.id, sim.name, sim.node);
			fscanf(mfptr, "%d%s%d", &sim.id, sim.name, &sim.node);
		}

		fclose(mfptr);
	}

	return 0;
}



