#ifndef UTILS
#define UTILS

#include "system.h"

void write_frame_diagnostics(std::ofstream& fout, double time, double Ekin, double Epot, double T_K);

namespace PhysicalConstants{
    constexpr double AVOGADRO {6.02214067e23} ;
    constexpr double K_B {8.314e-3 } ; // kJ/(mol * K)

}
#endif
