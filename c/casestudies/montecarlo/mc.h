
//System of spin 1/2 ising particles
#include <stdlib.h>
#include <time.h>

struct mc_particle
{

		int m_spin;

};


struct mc_system
{
		mc_particle *m_particles;	//Particle array
		int m_current;					//1d index of current particle
		int m_iside;					//length of square lattice
		int m_nparticles;				//number of particles
		int m_left;						//left neighbour to current particle
		int m_right;
		int m_up;
		int m_down;
		float m_h;

};


void init_mc_system(int iside, mc_system *m_system);
void perturb(mc_system *m_system);


int  find_index(int column, int row);
int  find_row(int i);
int  find_column(int i);

void find_neighbours(mc_system *m_system);

void sysenergy(mc_system *m_system);
void restore(mc_system *m_system);




