/*program23.c*/

/*Example of the use of malloc and calloc*/

/*malloc.c*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

struct siminfo
{
	int id;
	char name[30];
	int node;
};

typedef struct siminfo SIMINFO;
typedef SIMINFO * SIMINFOPTR;

main()
{
        int i; /*Control variable*/

	char nameroot[30] = "Default_sim_";
	char nname [30];
        int nid;
        int nnode;

	SIMINFOPTR simptr; /*Dynamically allocated  siminfo element*/

	SIMINFOPTR simarray;
	int isize; /*Array size variable*/

	srand(time(NULL)); /*Randomise random number generator*/

	printf("Using malloc to"
	       " dynamically allocate a single siminfo element.\n");
	simptr = (SIMINFOPTR)malloc(sizeof(SIMINFO));
	nid=1;
	nnode=1+rand()%50;
	
	simptr->id=nid;
	strcpy(simptr->name, "siminfoelement1");
	simptr->node=nnode;
	printf("\nThe single node is:\n");
	printf("%6s %20s %6s\n", "Runid", "Name", "Node");
	printf("%6d %20s %6d\n", simptr->id, simptr->name, simptr->node);

	printf("Enter a size for the array of structures.\n");
	scanf("%d", &isize);
	if(isize<1 || isize>500) isize=50;

	/*Create array of sim elements*/
	printf("Using calloc to create an array of %d sim nodes.\n", isize);
	simarray = (SIMINFOPTR)calloc(isize, sizeof(SIMINFO));
	for(i=0; i<isize; i++)
	{
		nid=i+1;
		nnode=1+rand()%50;
		sprintf(nname,"%s%d", nameroot, nid);

		simarray [i].id=nid;
		strcpy(simarray [i].name, nname);
		simarray [i].node=nnode;
	}

	/*Displaying the calloc'ed array*/
	printf("\nThe stored sim node array is.\n");
	printf("%6s %20s %6s\n", "Runid", "Name", "Node");
	for(i=0; i<isize; i++)
		printf("%6d %20s %6d\n", simarray[i].id, simarray[i].name, simarray[i].node);
	

        /*Avoid memory leaks!!!*/
        printf("Freeing the dynamically allocated memory.\n");
	free(simarray);
	free(simptr);

	return 0;



}
