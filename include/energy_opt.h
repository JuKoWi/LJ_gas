#ifndef ENERGY_OPT_H
#define ENERGY_OPT_H

#include "read_input.h"
#include "system.h"
#include <vector>

void energy_opt(std::vector<Particle>& system, Parameters& params);
void update_force_pbc(std::vector<Particle>& system, const Parameters& params);
std::vector<double> LJ_force(double dx, double dy, double dz, const Parameters& params);
void update_pos_pbc(std::vector<Particle>& system, const Parameters& params);
void velocity_verlet(std::vector<Particle>& system, const Parameters& params, double dt);
double pair_potential_energy(double dr2, Parameters params);
double potential_energy(const std::vector<Particle>& system, Parameters params);
#endif
