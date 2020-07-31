/*program16.c String Manipulation Library Functions*/

/*string.h */
/*strcpy
  strlen
  strcmp
  strtok
*/

/*stringfuncs.c*/

#include <stdio.h>
#include <string.h>


main()
{
	char sServer[] ="Bamford";
        char sAnotherServer[]="Derwent";
        char sCopy [20];
	char s1[] ="Titania,3.78 7"; /*sentence with 3 tokens*/

        char *tokenptr;
        char *tokenseparators = ", ";

        printf("Using the function strcpy(char *originalstring, char *copystring).\n");
	printf("The string in array sServer is: %s\n", sServer);
        printf("The copy of sServer sCopy is: %s\n\n", strcpy(sCopy, sServer));

        printf("\nUsing the strlen(char *s) function.\n");
	printf("sCopy %s is of length %d\n", sCopy, strlen(sCopy));
	printf("sServer %s is of length %d\n", sServer, strlen(sServer));

	printf("\nString Comparison Using strcmp(char *s1, char *s2)\n");
	printf("String sServer is %s\n", sServer);
	printf("String sAnotherServer is %s\n", sAnotherServer);
        printf("The copy string is %s\n", sCopy);

	printf("String sServer is the same as string sCopy:\n"); 
        printf("strcmp() returns %d\n", strcmp(sServer, sCopy)); 
	printf("String sServer is not the same as string sAnotherServer\n"); 
        printf("strcmp() returns %d\n", strcmp(sServer, sAnotherServer));	

	printf("\nUsing function strtok.\n");
        printf("The string to be tokenized is %s\n", s1);
	printf("The tokens are:\n");

	tokenptr = strtok(s1, tokenseparators);
	while(tokenptr != NULL)
	{
		printf("%s\n", tokenptr);
                tokenptr = strtok(NULL, tokenseparators);
	}

	return 0;
}
