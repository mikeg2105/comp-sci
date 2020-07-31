/*Program 2: Integer Input and Arithmetic*/
#include <stdio.h>

main()
{
	int i1, i2;                   /*Declare input variables*/
	int sum, difference;  /*Declare output variables*/

        /*Request input from the user*/
	printf("Enter the first integer\n"); 
	scanf("%d", &i1);       /*Read in the integer*/
	printf("Enter the second integer\n"); 
	scanf("%d", &i2); 

	/*Perform arithmetic operations*/ 
        sum=i1+i2;
        difference=i1-i2; 

	/*Display the results*/
        printf("The sum is %d\n", sum);
        printf("The difference between %d and %d is %d\n", i1, i2, difference);

	return 0;   /*Indicates the program ran successfully*/ 
}
