/*program14.c String input and output with sscanf and sprintf*/
/*stringio.c*/

#include <stdio.h>

main()
{
        char node[20], s2[80];
	char s1[] ="Titania 3.78 7";
        float fload, floadout;
        int nusers, nusersout;

        /*Using sscanf to read data from a string*/
        sscanf(s1, "%s%f%d", node, &floadout, &nusersout);
        puts("Using sscanf\n");
        printf("For node %s the statistics are:\n", node);
        printf("Load average=%f\n", floadout);
        printf("Number of users=%d\n", nusersout);
        
        /*Using sprintf to format a string*/
        puts("Enter node name");
        gets(node);
        printf("Enter load average and number of users for node %s:\n", node);
        scanf("%f%d", &fload, &nusers);

        sprintf(s2, "%s %f %d", node, fload, nusers);
        printf("The formatted output in s2 is:\n %s\n", s2);

	return 0;
}

