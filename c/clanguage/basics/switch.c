/*Program5 Using switch multiple-selection structures*/
/*switch.c*/

#include <stdio.h>

main()
{
	int file;
        int nMaxima=0, nTitania=0, nPascali=0, nSnowdon=0, nOther=0;

	printf("Enter the first letter of the node name (m, t, p or s).\n");
        printf("Enter the EOF character to end inut.\n");

        while((file=getchar()) != EOF)
        {
		/*switch nested in while loop*/

		/*This switch statement checks the upper an lower case
                  first character of the node name.*/
		switch(file)
		{
			case 'm': case 'M':
				++nMaxima;
				break;
			case 't': case 'T':
				++nTitania;
				break;
			case 'p': case 'P':
				++nPascali;
				break;
			case 's': case 'S':
				++nSnowdon;
				break;
			case '\n': case ' ':  /*Ignore these in input*/
				break;
			default:              /*Catch all other characters*/
				++nOther;
				break;
		} /*End of file check switch */
	} /*End of file check while loop*/

	printf("\nThe totals for each node are:\n");
	printf("Titania: %d\n", nTitania);
	printf("Maxima: %d\n", nMaxima);
	printf("Pascali: %d\n", nPascali);
	printf("Snowdon: %d\n", nSnowdon);
	printf("Other: %d\n", nOther);

}
