#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main()
{
	int row,column;

	srand((unsigned)time(NULL));

	for(row=0;row<10;row++)
	{
		for(column=0;column<10;column++)
			printf("%2d\t",rand() % 100);
		putchar('\n');
	}

	return(0);
}

