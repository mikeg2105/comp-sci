#include <stdio.h>
#include <string.h>

int main()
{
	char password[9];
	int compare;

	printf("Enter your password: ");
	scanf("%8s",password);
	compare = strcmp(password,"secret");
	if( compare == 0 )
		puts("You've been granted access!");
	else
		puts("Intruder alert!");

	return(0);
}

