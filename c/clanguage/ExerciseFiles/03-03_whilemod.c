#include <stdio.h>

int main()
{
	int a,b;

	a = 7;
	while( a < 30 )
	{
		b = a % 7;		/* b equals a mod 7 */
		printf("%d %% 7 = %d\n",a,b);
		a++;
	}

	return(0);
}

