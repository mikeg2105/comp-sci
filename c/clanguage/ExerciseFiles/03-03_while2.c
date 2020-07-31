#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main()
{
	int x,r;

	srand((unsigned)time(NULL));	/* Seed randomizer */
	
	x = 10;
	while( x > 0 )
	{
		r = rand();
		printf("%d\n",r);
		x--;
	}

	return(0);
}

