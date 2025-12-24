#include <vector>
#include <cmath>
#include <string>
#include "read_input.h"

struct Particle{
	double x, y, z;
	double vx, vy, vz;
	double fx, fy, fz;
};

void initialize_system(std::vector<Particle>& system, Parameters params){
	for (auto& p : system){
		p.x = drand48() *  params.pbcx_angstrom - 0.5 * params.pbcx_angstrom;
		p.y = drand48() *  params.pbcy_angstrom - 0.5 * params.pbcy_angstrom;
		p.z = drand48() *  params.pbcz_angstrom - 0.5 * params.pbcz_angstrom;
		p.vx = drand48() * 9.0 -4.5;
		p.vy = drand48() * 9.0 -4.5;
		p.vz = drand48() * 9.0 -4.5;
	}
}
