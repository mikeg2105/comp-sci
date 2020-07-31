/*program19.c The Enumeration type*/
/*enum.c*/

#include <stdio.h>

enum bool {true=1, false=0, TRUE=1, FALSE=0, yes=1, no=0, YES=1, NO=0, 
           ON=1, OFF=0, on=1, off=0};

enum months {JAN=1, FEB, MAR,APR,MAY,JUN,JUL,AUG,SEP,OCT,NOV,DEC};

typedef enum bool boolean;


boolean testfunc();

main()
{
        enum months month;

        char *sMonth[] ={"", "January", "February", "March", "April", 
                             "May", "June", "July", "August", "September",
                             "October", "November", "December"};

        for(month=JAN; month<=DEC; month++)

	        /*The format specifier used in printf
                  specifies the width of the integer output
                  and the width of the string output */
		printf("%2d%11s\n", month, sMonth[month]);


	printf("You hit %s!\n", (booleantestfunc()?"yes":"no"));

	return 0;
}

booleantestfunc()
{
	boolean btest=false;
        boolean finished=false;
	char c;

	while(!finished)
	{
		printf("Enter Yes(y) or No(n).\n");
                finished=true;
		switch(c=getchar())
		{
			case 'y': case 'Y':
			  	btest=true;
			break;
			case 'n': case 'N':
				btest=false;
			break;
			default:
				finished=false;
			break;
		}		
	}

	return btest;
}
