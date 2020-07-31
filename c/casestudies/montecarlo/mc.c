double randno()
{
	return 1.0*rand()/(pow(2,31)-1);		
}


init_mc_system(int iside, mc_system *m_system)
{
                m_system->m_iside=iside;		
		m_system->m_nparticles=iside*iside;
		m_system->m_particles=   (mc_particle *)calloc(m_system->m_nparticles,sizeof(mc_particle));
                for(int i=0; i<m_system->m_nparticles; i++)
                     m_system->m_particles[i].m_spin=1.0;
		m_system->m_h=0.1;							
}

perturb(mc_system *m_system)
{
		m_system->m_current=randno();
		flipspin(m_system);
}

flipspin()
{		
			m_system->m_particles[m_system->m_current]=-m_system->m_particles[m_system->m_current];
}



int  find_index(int column, int row){return m_iside*row+column;}
int  find_row(int i){return i/m_iside;}
int  find_column(int i){return i%m_iside;}

void find_neighbours(mc_system *m_system)
{
	int left_column, right_column, up_row, down_row;
	
	int column=find_column(m_system->m_current);
	int row=find_row(m_system->m_current);
	
	if(column==0)
		left_column=(m_system->m_iside) -1;
	else
		left_column=column-1;
		
	if(column==(m_system->m_iside))
		right_column=0;
	else
		right_column=column+1;
		
	if(row==0)
		up_row=(m_system->m_iside)-1;
	else
		up_row=row-1;
		
	if(row==(m_system->m_iside)-1)
		down_row=0;
	else
		down_row=row+1;
		
	m_system->m_left=find_index(left_column, row);
	m_system->m_right=find_index(right_column, row);
	m_system->m_up=find_index(column, up_row);
	m_system->m_down=find_index(column,down_row);
		
}

double sysenergy(mc_system *m_system)
{
		find_neighbours();
                mc_particle *particles=m_system->m_particles;

                mc_particle mcurrent=particles[m_current];
		mc_particle mleft=particles[m_system->m_left];
                mc_particle mright=particles[m_system->m_right];
                mc_particle mup=particles[m_system->m_up];
                mc_particle mdown=particles[m_system->m_down];
      
		
		double E= -(mcurrent.m_spin)*(    mleft.m_spin+mright.m_spin     );

		E-=mcurrent*(mup+mdown);
		
		return 2.0*E;
}

void restore(mc_system *m_system)
{
			m_system->m_particles[m_system->m_current]=-m_system->m_particles[m_system->m_current];
}

double magnetisation(mc_system *m_system)
{
	double mag=0;
        mc_particles particles=m_system.m_particles;
	
	for(int i=0; i<(m_system->m_nparticles); i++)
        {
		mag+=particles[i].m_spin;
        }
	return mag/(m_system->m_nparticles);	
}



