#include <stdio.h>

int main()
{
	int pokey;
	int *p;

	pokey = 987;
	p = &pokey;
	
	printf("The address of `pokey` is %p\n",&pokey);
	printf("The contents of `pokey` are %d\n",pokey);

	printf("The address of `pokey` is %p\n",p);
	printf("The contents of `pokey` are %d\n",*p);

	return(0);
}

