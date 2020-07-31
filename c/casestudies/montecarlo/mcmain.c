//

#include "mc.h"


#include <stdio.h>



int main(int argc, char **argv)
{
	int Nside, Niterations;
	float temp;

        FILE *fout;

        if((fout=fopen("results.txt"),"w")==NULL)
        {
            printf("File results.txt not opened\n");
            return 1;
        }
	
	printf( "This is an Ising Particle Monte Carlo Simulation\n" );
	
	scanf( "Enter the length of the lattice: %d\n",&Nside);	
	scanf( "Enter the number of Monte-carlo iterations: ", &Niterations);	
	//scanf( "Enter the temperature: \n", &temp);
	
	/* initialize random seed: */
	srand ( time(NULL) );
	
	for(int i=1; i<30; i++)
	{
		temp = 0.1*i;
		MC_simulation(Nside, Niterations, temp, fout);
	}
	
	fclose(fout);
	return 0;	
}

void MC_simulation(int Nside, int Niterations, double T, FILE *fout)
{
	//Instantiate a two dimensional mc system of square Nside
        mc_system mcsystem;
        float tmagnetisation=0;
        float tmagsq=0;
        int ttotal=0;
        float variance;

        init_mcsystem(Nside, &mc_system);
	
	
	
	//Warm up iterations
	for(int iter=1; iter<=100*Nside*Nside; iter++)
	{
		perturb(&mc_system);
		double energy=sysenergy(&mc_system);
		if(exp( -energy/T)<(rand()%100)/100)
					restore(&mc_system);		
	}
	
	//The Monte-carlo simulation loop
	for(int iter=1; iter<=Niterations; iter++)
	{
		perturb(&mc_system);
		double energy=sysenergy(&mc_system);
		if(exp( -energy/T)<(randno())
					restore(&mc_system);
		else
		{
			double magnetisation=magnetisation(&mc_system);
			tmagnetisation +=magnetisation;
                        tmagsq +=magnetisation*magnetisation
                        ttotal++;	
		}			
	}
	
        variance=(ttotal*tmagsq-tmagnetisation*tmagnetisation)/((ttotal-1)*ttotal);
        variance=sqrt(variance/ttotal);

	fprintf(fout,"Monte Carlo simulation of a  %d x %d  spin 1/2 mc system.\n",Nside,Nside);
	fprintf(fout,"The temperature is: %f\n" , T );
	fprintf(fout,"The number of iterations were: %d\n", Niterations );
	fprintf(fout,"The average magnetisation is: %f\n", tmagnetisation/ttotal);
	fprintf(fout,"The standard error of the mean is: %f\n", variance);
}
