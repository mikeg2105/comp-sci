#include <stdio.h>
#include <stdlib.h>

struct stats {
	int account;
	float balance;
	struct stats *next;
};

struct stats *create(void);

int main()
{
	struct stats *first;
	struct stats *current;
	struct stats *new;		/* Change name if necessary */
	int a = 1;

	/* create first structure */
	first = create();

	return(0);
}

/* Allocate storage for one new structure */
struct stats *create(void)
{
	struct stats *baby;

	baby = (struct stats *)malloc(sizeof(struct stats));
	if( baby == NULL)
	{
		puts("Memory error");
		exit(1);
	}
	return(baby);
}

