#include <stdio.h>

int main()
{
	char *gang[5] = {
		"Spanky",
		"Buckwheat",
		"Alfalfa",
		"Darla",
		"Pete the Pup"
	};
	int x;
	char *cptr;

	for(x=0;x<5;x++)
	{
		cptr = gang[x];
		while(*cptr)
		{
			putchar(*cptr);
			cptr++;
		}
		putchar('\n');
	}

	return(0);
}

