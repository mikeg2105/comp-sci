/*program18.c*/
/* Illustration of typedef and using structures with functions*/

/*typedef.c*/

#include <stdio.h>

struct node {
	char *name;
	char *processor;
        int  num_procs;
	};


typedef struct node Node;

void FillGrid(Node *, char *[], char *[], int []);
void DisplayGrid(Node *);


main()
{
   Node wrgrid[4];

   char *node_names[]={"Titania", "Maxima", "Pascali", "Snowdon"};
   char *proc_names[]={"Ultra Sparc III Cu","Intel Pentium Xeon"};
   int  num_procs[]={80, 40, 8, 256};

   FillGrid(wrgrid, node_names, proc_names, num_procs);
   DisplayGrid(wrgrid);

   return 0;
}

/*Fill a grid definition with information
  about nodes*/
void FillGrid(Node *Grid, char *sNames[], char *sProcs[], int nProcs[])
{
	int i;
	for(i=0; i<4; i++)
	{
		Grid[i].name=sNames [i];
                Grid[i].processor=(i==3?sProcs[1]:sProcs[0]); /*Note use of  */
                                                              /* trinary (?:)*/
                                                              /*  operator   */
	        Grid[i].num_procs=nProcs[i];
	}
}

/*Display the contents for each node*/
void DisplayGrid(Node *Grid)
{
	int i;
	printf("The White Rose Grid Nodes are:\n");
	for(i=0;i<4;i++)
           printf("%s with %d %s processors.\n",
		 Grid[i].name, Grid[i].num_procs, Grid[i].processor);

}

