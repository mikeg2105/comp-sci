/*Example4.c:Counter-controlled repetition with the for structure*/
/*for.c*/
#include <stdio.h>

main()
{
	int counter, nsteps=10;
        
	/* initialisation, repetition condition and increment */
        /* are all included in the for structure header       */
        for(counter=1; counter<=nsteps; counter++)
			printf("%d\n", counter);

	return 0;

}
