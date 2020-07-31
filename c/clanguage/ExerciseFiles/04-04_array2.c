#include <stdio.h>

int main()
{
	float temps[4] = { 84.9, 83.7, 85.8, 88.2 };
	int x;

	printf("Local temperatures:\n");
	for(x=0;x<4;x++)
		printf("Station %d: %.1f\n",x,temps[x]);

	return(0);
}

