/*program24.c*/
/*Linked list example*/

/*linkedlist.c*/

#include <stdio.h>
#include <stdlib.h>

struct siminfo  /*Self referential structure*/
{
	int id;
	char name[30];
	int node;
	struct siminfo *nextptr;
};

typedef struct siminfo SIMINFO;
typedef SIMINFO *SIMINFOPTR;

void insert(SIMINFOPTR *, int, char *, int);
int delete(SIMINFOPTR *, int);
int isempty(SIMINFOPTR );
void displaylist(SIMINFOPTR );
void displayhelp();

main()
{
	SIMINFOPTR startptr = NULL;
        int choice;

        int in_id;
        char in_name[30];
        int in_node;

	displayhelp(); /*Display help menu*/
	printf("enter choice? ");
	scanf("%d", &choice);

	while(choice != 5)
	{
		switch(choice)
		{
			case 1: 
				displayhelp();
				break;
			case 2:
				printf("Enter the runid, simulation name, and node id.\n");
				scanf("%d%s%d", &in_id, in_name, &in_node);
				insert(&startptr, in_id, in_name, in_node);
				displaylist(startptr);
				break;
		
			case 3:
				if(!isempty(startptr))
				{
					printf("Enter ID to be deleted:");
					scanf("\n%d", &in_id);

					if(delete(&startptr, in_id))
					{
						printf("node id %d deleted.\n", in_id);
						displaylist(startptr);
					}
					else
						printf("node id %d not found.\n", in_id);
				}
				else
					printf("The node list is empty.\n");

				break;
			case 4:
				displaylist(startptr);
				break;
			default:
				printf("Invalid choice.\n");
				displayhelp();
				break;
		}
		printf("enter choice? ");
		scanf("%d", &choice);
	}

	printf("End of simnode linked list run\n");
	return 0;
}

/*insert value into the list in order of id*/
void insert(SIMINFOPTR *sptr, int vid, char *vname, int vnode)
{
	SIMINFOPTR newptr, previousptr, currentptr;

	newptr = malloc(sizeof(SIMINFO));
	
	if(newptr != NULL) /*space is available*/
	{
		newptr->id=vid;
		strcpy(newptr->name, vname);
		newptr->node=vnode;
		newptr->nextptr=NULL;

		previousptr=NULL;
		currentptr=*sptr;

                /*Walk to correct node in linked list*/
		while(currentptr!=NULL && vid > currentptr->id)
		{
			previousptr=currentptr;
			currentptr=currentptr->nextptr;
		}

		if(previousptr==NULL)
		{
			newptr->nextptr=*sptr;
			*sptr=newptr;
		}
		else
		{
			previousptr->nextptr=newptr;
			newptr->nextptr=currentptr;
		}
	}
	else
	    printf("%d %s %d not inserted, out of memory!\n", vid, vname, vnode);
	

}


/*Delete a list element*/
int delete(SIMINFOPTR *sptr, int vid)
{
	SIMINFOPTR previousptr, currentptr, tempptr;

	if(vid==(*sptr)->id)
	{
		tempptr=*sptr;
		*sptr=(*sptr)->nextptr; /*remove node from list*/
		free(tempptr);          /*free the removed node*/
		return vid;
	}
	else
	{
		previousptr= *sptr;
		currentptr= (*sptr)->nextptr;

		/*Walk to node to be deleted*/
		while(currentptr!=NULL && currentptr->id != vid)
		{
			previousptr=currentptr;
			currentptr=currentptr->nextptr;
		}

		if(currentptr != NULL)
		{
			tempptr=currentptr;
			previousptr->nextptr=currentptr->nextptr;
			free(tempptr);
			return vid;
		}
	}

	return 0;
}

/*if the list is empty return 1*/
int isempty(SIMINFOPTR sptr )
{
	return sptr==NULL;
}

void displaylist(SIMINFOPTR currentptr)
{
	if(currentptr==NULL)
		printf("The list is empty\n");
	else
	{
		printf("The list nodes are:\n");	

		while(currentptr != NULL)
		{
			printf("%6d %20s %6d\n", currentptr->id, currentptr->name, currentptr->node);
			currentptr=currentptr->nextptr;
		}
		printf("NULL\n\n");
	}
}

void displayhelp()
{
	printf("Enter a menu selection:\n"
		"	1 to display help.\n"
		"	2 to insert an element into the list.\n"
		"	3 to delete an element from the list.\n"
		"	4 to display the node list.\n"
		"	5 to end.\n");
}

