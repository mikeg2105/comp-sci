#include <stdio.h>

int main()
{
	char *string = "I'm just a normal string.\n";
	int x = 0;

	while(string[x])
	{
		putchar(string[x]);
		x++;
	}

	return(0);
}

