#ifndef ENERGY_OPT_H
#define ENERGY_OPT_H

#include "read_input.h"
#include "system.h"
#include <vector>

void energy_opt(System& system, const Parameters& params);
void update_force_pbc(System& system, const Parameters& params, bool verlet = true);
std::vector<double> LJ_force(double dx, double dy, double dz, const Parameters& params);
void update_pos_pbc(System& system);
void velocity_verlet(System& system, const Parameters& params, double dt);
double pair_potential_energy(double dr2, const Parameters params);
double potential_energy(const System& system, const Parameters params);
double kinetic_energy(const System& system);
double minimum_image_distance(double dx, double L);
#endif
