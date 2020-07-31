#ifndef _MC_PARTICLE_H_
#define _MC_PARTICLE_H_


class Cmc_Particle
{
	private:
		int m_spin;
	
	
	public:
		Cmc_Particle();
		~Cmc_Particle();
		void flip_spin(){m_spin=-m_spin;}
		inline int spin_value(){return m_spin;}	
	
	
	
};

#endif //_mc_PARTICLE_H_
