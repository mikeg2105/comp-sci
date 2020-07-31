#include <stdio.h>

int main()
{
	char hello[] = "Greetings, human!\n";
	int n;

	n = 0;
	while( putchar(hello[n++]) )	/* ugly, huh? But it works! */
		;							/* This is common in C when an
									  empty statement belongs to
   									  a loop or `if` decision. It
									  shows that the empty statement
									  is intentional.
									 */

	return(0);
}

