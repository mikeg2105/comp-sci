#ifndef _MC_STATISTICS_H_
#define _MC_STATISTICS_H_
#include <iostream>
#include <cmath>
#include <cstdlib>

double randno();

class CStatistics
{
	private:
		double m_sum_x;
		double m_sum_xsq;
		long int m_N;
	
	public:
		CStatistics();
		void add(double x);
		double average();
		double standard_error();	
};



#endif //_MC_STATISTICS_H_
