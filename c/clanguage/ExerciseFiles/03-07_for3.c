#include <stdio.h>

int main()
{
	char ch,t;

	printf("Type a letter 'a' to 'z': ");
	scanf("%c",&ch);

	for(t='a';t<='z';t++)
	{
		putchar(t);
	}
	putchar('\n');

	return(0);
}

