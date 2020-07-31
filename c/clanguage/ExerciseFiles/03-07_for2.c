#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main()
{
	int column,r;

	srand((unsigned)time(NULL));	/* Seed randomizer */

	for(column=0;column<10;column++)
	{
		r = rand();
		printf("%2d\t",r % 100);    /* keep random numbers between 0 and 99 */
	}
	putchar('\n');                  /* newline output, end of line */

	return(0);
}

