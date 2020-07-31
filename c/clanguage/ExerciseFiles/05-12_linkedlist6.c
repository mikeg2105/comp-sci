#include <stdio.h>
#include <stdlib.h>

struct stats {
	int account;
	float balance;
	struct stats *next;
};

void fill_structure(struct stats *s);
struct stats *create(void);

int main()
{
	struct stats *first;
	struct stats *current;
	struct stats *new;		/* Change name if necessary */
	int x = 5;

	/* create first structure */
	first = create();
	current = first;

	for(x=0;x<5;x++)
	{
	    if(x==0)
        {

            first = create();
            current = first;
        }
        else
        {
            new = create();
            current->next = new;
            current = new;
        }
		fill_structure(current);
	}
	current->next = NULL;

	current = first;		/* reset the list */
	while(current)
	{
		printf("Account %05d:\t$%.2f\n",
				current->account,
				current->balance);
		current = current->next;
	}

	return(0);
}

/* Fill a structure */
void fill_structure(struct stats *s)
{
	printf("Account number: ");
	scanf("%d",&s->account);
	printf("Balance: $");
	scanf("%f",&s->balance);
	s->next = NULL;
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

