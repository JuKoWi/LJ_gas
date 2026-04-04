#ifndef UTILS
#define UTILS

#include "system.h"

void write_xyz(std::ofstream& fout, const std::vector<Particle>& system);
void write_energy_force(std::ofstream& fout, double time, double Ekin, double Epot);

namespace PhysicalConstants{
    constexpr double AVOGADRO = 6.02214067e23;
}

#endif
