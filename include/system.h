#ifndef SYSTEM_H
#define SYSTEM_H

#include "read_input.h"
#include <vector>

struct Particle{
	double x, y, z; // nm
	double vx, vy, vz; // nm/ps
	double fx, fy, fz;
};

void initialize_random_positions(std::vector<Particle>& system, Parameters params);
void initialize_positions_sc(std::vector<Particle>& system, Parameters params);
void initialize_constant_velocity(std::vector<Particle>& system);
void initialize_rand_velocity(std::vector<Particle>& system);
void initialize_positions_file(std::vector<Particle>&);
#endif
