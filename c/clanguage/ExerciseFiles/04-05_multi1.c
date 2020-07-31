#include <stdio.h>

#define ROWS 4
#define COLUMNS 4

int main()
{
	int grid[ROWS][COLUMNS];
	int x,y;

	/* initialize the array */
	for(x=0;x<ROWS;x++)
		for(y=0;y<COLUMNS;y++)
			grid[x][y] = 0;

	/* display the grid */
	for(x=0;x<ROWS;x++)
	{
		for(y=0;y<COLUMNS;y++)
			printf("%d.%d: %d\t",x,y,grid[x][y]);
		putchar('\n');
	}

	return(0);
}

