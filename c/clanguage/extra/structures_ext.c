/*program17.c */

/*structures.c*/

/*Using and defining structures*/
#include <string.h>

struct node {
	char *name;
	char *processor;
        int  num_procs;
	};

typedef struct node NODE;
typedef struct node *NODEPTR;

NODEPTR addnode(NODE node1, NODE node2);
main()
{
	struct node n1;
	struct node n2;
	struct node *n1ptr;
	struct node *addptr;
	
	
	n1.name="Titania";
	n1.processor ="Ultra Sparc III Cu";
	n1.num_procs = 80;
	
	n2.name="Iceberg";
	n2.processor ="AMD opteron";
	n2.num_procs = 900;
	
	n1ptr = &n1;

	printf("Access using the structure\n");
	printf("The node %s has %d %s processors\n", 
                n1.name, n1.num_procs, n1.processor);

	printf("Access using the structure pointer\n");
	printf("The node %s has %d %s processors\n", 
                n1ptr->name, n1ptr->num_procs, n1ptr->processor);

	printf("Access by dereferncing the structure pointer\n");
	printf("The node %s has %d %s processors\n", 
                (*n1ptr).name, (*n1ptr).num_procs, (*n1ptr).processor);

    addptr=addnode(n1,n2);
 	printf("Result of adding two structures\n");
	printf("The addnode %s has %d %s processors\n", 
                (*addptr).name, (*addptr).num_procs, (*addptr).processor);
   
    free(addptr);
    //system("PAUSE");
	return 0;
}

NODEPTR addnode(NODE node1, NODE node2)
{
       struct node *n1ptr;       
       n1ptr=(NODEPTR)malloc(sizeof(NODE));  
       n1ptr->name="added processors";
       n1ptr->processor="compound";  
       n1ptr->num_procs=(node1.num_procs)+(node2.num_procs);
       return n1ptr;
}


