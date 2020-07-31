/*program17.c */

/*structures.c*/

/*Using and defining structures*/

struct node {
	char *name;
	char *processor;
        int  num_procs;
	};

main()
{
	struct node n1;
	struct node *n1ptr;

	n1.name="Titania";
	n1.processor ="Ultra Sparc III Cu";
	n1.num_procs = 80;

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


	return 0;
}
