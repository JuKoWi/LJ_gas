#ifndef SYSTEM_H
#define SYSTEM_H

#include "read_input.h"
#include <vector>

struct Particle{
	double x, y, z;
	double vx, vy, vz;
	double fx, fy, fz;
};

void initialize_random_positions(std::vector<Particle>& system, Parameters params);
void initialize_positions_sc(std::vector<Particle>& system, Parameters params);
#endif
