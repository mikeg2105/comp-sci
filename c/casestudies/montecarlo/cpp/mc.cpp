
//1st April 2005
//Taken from a C++ teaching example by:
//William Barford and Lee Thomson
//Department of Physics
//The University of Sheffield

#include "../include/mc.hpp"
#include "../include/mc_statistics.hpp"

using std::fstream;
using std::cout;
using std::endl;
using std::cin;

int main(int argc, char **argv)
{
	fstream fout;
	
	int Nside, Niterations;
	float temp;
	fout.open("results.txt", std::ios::out);
	
	cout << "This is an Ising Particle Monte Carlo Simulation" << endl;
	
	cout << "Enter the length of the lattice: ";
	cin >> Nside;
	
	cout << "Enter the number of Monte-carlo iterations: ";
	cin >> Niterations;
	
	//cout << "Enter the temperature: ";
	//cin >> temp;
	
	for(int i=1; i<30; i++)
	{
		temp = 0.1*i;
		MC_simulation(Nside, Niterations, temp, fout);
	}
	
	fout.close();
	return 0;	
}

void MC_simulation(int Nside, int Niterations, double T, fstream &fout)
{
	//Instantiate a two dimensional mc system of square Nside
	Cmc_System MC_system(Nside);
	CStatistics MC_statistics;
	
	//Warm up iterations
	for(int iter=1; iter<=100*Nside*Nside; iter++)
	{
		MC_system.perturb();
		double energy=MC_system.energy();
		if(exp( -energy/T)<randno())
					MC_system.restore();		
	}
	
	//The Monte-carlo simulation loop
	for(int iter=1; iter<=Niterations; iter++)
	{
		MC_system.perturb();
		double energy=MC_system.energy();
		if(exp( -energy/T)<randno())
					MC_system.restore();
		else
		{
			double magnetisation=MC_system.magnetisation();
			MC_statistics.add(magnetisation);	
		}			
	}
	
	fout<<"Monte Carlo simulation of a " << Nside << "x" << Nside << " spin 1/2 mc system." << endl;
	fout<<"The temperature is: " << T << endl;
	fout<<"The number of iterations were: " << Niterations << endl;
	fout<<"The average magnetisation is: " << MC_statistics.average()<<endl;
	fout<<"The standard error of the mean is: " << MC_statistics.standard_error()<<endl<<endl;
}
