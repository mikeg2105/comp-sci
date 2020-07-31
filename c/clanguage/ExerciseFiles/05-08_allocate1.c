#include <stdio.h>

int main()
{
	char c;
	int i;
	float f;
	char buffer[128];
	long unsigned total;

	total = 0;
	total = total + sizeof(c);
	total = total + sizeof(i);
	total = total + sizeof(f);
	total = total + sizeof(buffer);
	total = total + sizeof(total);

	printf("This program set aside %ld bytes of storage.\n",total);

	return(0);
}

