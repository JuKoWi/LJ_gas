#ifndef SYSTEM_H
#define SYSTEM_H

struct Particle{
	double x, y, z;
	double vx, vy, vz;
	double fx, fy, fz;
};

void initialize_system(std::vector<Particle>& system, Parameters params);

#endif
