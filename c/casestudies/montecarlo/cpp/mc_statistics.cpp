	#include "../include/mc_statistics.hpp"

	double ::randno()
	{
		return 1.0*rand()/(pow(2,31)-1);		
	}
	
	CStatistics::CStatistics()
	{
		m_sum_x=0.0;
		m_sum_xsq=0.0;
		m_N=0;
	}
	
	void CStatistics::add(double x)
	{
		m_sum_x+=x;
		m_sum_xsq+=x*x;
		m_N++;	
	}
	
	double CStatistics::average()
	{
		return m_sum_x/m_N;
	}
	
	double CStatistics::standard_error()
	{
		double variance=(m_N*m_sum_xsq-m_sum_x*m_sum_x)/((m_N-1)*m_N);
		return sqrt(variance/m_N);	
	}
	
