#ifndef _MC_H_
#define _MC_H_

#include <iostream>
#include <fstream>
#include "mc_system.hpp"

using std::fstream;

void MC_simulation(int Nside, int Niterations, double T, fstream &fout);

#endif //_MC_H_
