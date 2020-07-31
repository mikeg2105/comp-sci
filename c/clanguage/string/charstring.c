/*program13.c String definition and input using gets, getchar and puts*/

/*charstring.c*/
#include <stdio.h>

main()
{
	char node[] = "Titania";
        char *nodeptr = "Titania";

        /*Alternative array definition of a string*/
        /*IMPORTANT NOTE the terminating NULL character.*/
	/*char node[] = {'T', 'i', 't', 'a', 'n', 'i', 'a', '\0'};*/

        char c, nodename1[80], nodename2[80];
        int i=0;

	puts("Enter a line of text");
	while((c=getchar())!='\n')
			nodename1[i++]=c;
	nodename1[i]='\0';
        printf("The line entered was %s\n", nodename1);

	printf("Enter another line of text\n");
        gets(nodename2);
	puts("The line entered was: ");
        puts(nodename2);

	return 0;

}
