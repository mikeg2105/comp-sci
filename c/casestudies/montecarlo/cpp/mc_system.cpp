
	#include "../include/mc_system.hpp"

	Cmc_System::Cmc_System(int iside)
	{
		m_iside=iside;
		m_Nparticles=m_iside*m_iside;
		m_particles=new Cmc_Particle [m_Nparticles];
		m_h=0.1;							
	}
		
	Cmc_System::~Cmc_System()
	{
		delete [] m_particles;
		cout << "Monte-carlo simulation completed. Good bye." << endl<<endl;			
	}
	
	void Cmc_System::perturb()
	{
		m_current=m_Nparticles*randno();
		m_particles[m_current].flip_spin();
	}
	
	void Cmc_System::restore()
	{
		m_particles[m_current].flip_spin();	
	}
	
	double Cmc_System::energy()
	{
		find_neighbours();
		
		double E= -m_particles[m_current].spin_value()*(m_particles[m_left].spin_value()+m_particles[m_right].spin_value());
		E-=m_particles[m_current].spin_value()*(m_particles[m_up].spin_value()+m_particles[m_down].spin_value());
		
		return 2.0*E;
	}
	
	double Cmc_System::magnetisation()
	{
		double mag=0;
		
		for(int i=0; i<m_Nparticles; i++)
			mag+=m_particles[i].spin_value();
		return mag/m_Nparticles;	
	}
	
	void Cmc_System::find_neighbours()
	{
		int left_column, right_column, up_row, down_row;
		
		int column=find_column(m_current);
		int row=find_row(m_current);
		
		if(column==0)
			left_column=m_iside -1;
		else
			left_column=column-1;
			
		if(column==m_iside)
			right_column=0;
		else
			right_column=column+1;
			
		if(row==0)
			up_row=m_iside-1;
		else
			up_row=row-1;
			
		if(row==m_iside-1)
			down_row=0;
		else
			down_row=row+1;
			
		m_left=find_index(left_column, row);
		m_right=find_index(right_column, row);
		m_up=find_index(column, up_row);
		m_down=find_index(column,down_row);
			
	}
	
	


	
	
	

		