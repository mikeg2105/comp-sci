#include <stdio.h>
#include <ctype.h>

int main()
{
	int ch = 'a';

	printf("Big %c\n",toupper(ch));
	printf("Little %c\n",tolower(ch));
	printf("What begins with %c?\n",ch);

	return(0);
}

