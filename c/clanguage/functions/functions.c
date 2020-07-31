/*Program7  Function calls*/
/*functions.c*/
/*A simple dice game*/
/*Each player starts with a random number of points*/
/*Each player takes turn to roll the dice
  the score is removed form the opposing players
  points*/

/*The game is won when a player has zero or less points*/
/*The winning player has more greater than zero points*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int RollDice(int iPlayer);
int SetPlayerPoints(int iPlayer);

main()
{
	int gamestatus=1;
        int p1points, p2points;

	/*Set player points*/
	p1points=SetPlayerPoints(1);
        p2points=SetPlayerPoints(2);

	while(gamestatus==1)
	{
		if(p1points>0) p2points -= RollDice(1);
		if(p2points>0) p1points -= RollDice(2);

		if((p1points<=0) || (p2points<=0)) gamestatus=0;
		printf("\n Player 1 has %d points.\n", p1points);
		printf("Player 2 has %d points.\n", p2points);
	}

	if(p1points>p2points)
		printf("Player1 wins!!\n");
	else
		printf("Player2 wins!!\n");

	return 0;
}


/*Roll a pair of 6 sided dice for player with
  id given by i player return the total */

int RollDice(int iPlayer)
{
	int d1, d2;
	int sum;

	d1 = 1+(rand()%6);
	d2 = 1+(rand()%6);
	sum = d1+d2;

	if((d1==d2) && d1==1) sum*=-1;
	if((d1==d2) && d1==6) sum*=2;

	printf("Player %d rolled %d and %d = %d\n", iPlayer, d1, d2, sum);

	return sum;
}

int SetPlayerPoints(int iPlayer)
{
	int iPoints;

        iPoints = 1+(rand()%100);
	printf("Player %d has %d points\n", iPlayer, iPoints);

	return iPoints;
}
