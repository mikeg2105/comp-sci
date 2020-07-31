#include <stdio.h>
#include <stdlib.h>

int main()
{
	struct date {
		int month;
		int day;
		int year;
	};
	struct date *birthday;

	birthday = (struct date *)malloc(sizeof(struct date));
	if( birthday == NULL)
	{
		puts("Memory unavailable");
		return(1);
	}

	printf("Enter your birth year: ");
	scanf("%d",&birthday->year);
	printf("Enter your birth month: ");
	scanf("%d",&birthday->month);
	printf("Enter your birth day: ");
	scanf("%d",&birthday->day);

	printf("You were born on %d/%d/%d\n",
			birthday->month,
			birthday->day,
			birthday->year);

	return(0);
}

