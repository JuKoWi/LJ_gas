#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <string>
#include "read_input.h"
#include "system.h"
#include "utils.h"
//comment


void apply_pbc(std::vector<Particle>& system, const Parameters& params){ // 0-centered simulation box
	for (auto& p : system){
		if (p.x > params.pbcx_angstrom*0.5){
			p.x -= params.pbcx_angstrom * std::floor(p.x / params.pbcx_angstrom);
		}
		if (p.y > params.pbcy_angstrom*0.5){
			p.y -= params.pbcy_angstrom * std::floor(p.y / params.pbcy_angstrom);
		}
		if (p.z > params.pbcz_angstrom*0.5){
			p.z -= params.pbcz_angstrom * std::floor(p.z / params.pbcz_angstrom);
		}
		if (p.x < -params.pbcx_angstrom*0.5){
			p.x -= params.pbcx_angstrom * std::floor(p.x / params.pbcx_angstrom);
		}
		if (p.y < -params.pbcy_angstrom*0.5){
			p.y -= params.pbcy_angstrom * std::floor(p.y / params.pbcy_angstrom);
		}
		if (p.z < -params.pbcz_angstrom*0.5){
			p.z -= params.pbcz_angstrom * std::floor(p.z / params.pbcz_angstrom);
		}
	}
}

std::vector<double> LJ_force(double dx, double dy, double dz, const Parameters& params)
{
    double r2 = dx*dx + dy*dy + dz*dz;
    if (r2 < 1e-12) return {0.0, 0.0, 0.0};

    double inv_r2 = 1.0 / r2;
    double inv_r6 = inv_r2 * inv_r2 * inv_r2;
    double inv_r12 = inv_r6 * inv_r6;

    double f = 24 * params.epsilon * inv_r2 * (2 * std::pow(params.sigma,12) * inv_r12 - std::pow(params.sigma, 6) * inv_r6);

    return {f * dx, f * dy, f * dz};
}

void update_force_pbc(std::vector<Particle>& system, const Parameters& params){
	for (auto& p : system){
		p.fx = 0;
		p.fy = 0;
		p.fz = 0;
	}

	for (int i = 0; i<params.N; ++i){
		for (int j = i + 1; j < params.N; ++j){ // only iterate over pairs
			double dx = system[j].x - system[i].x;
			double dy = system[j].y - system[i].y;
			double dz = system[j].z - system[i].z;
			if (dx > 0.5 * params.pbcx_angstrom){
				dx -= params.pbcx_angstrom;
			}
			if (dx < -0.5 * params.pbcx_angstrom){
				dx += params.pbcx_angstrom;
			}
			if (dy > 0.5 * params.pbcy_angstrom){
				dy -= params.pbcy_angstrom;
			}
			if (dy < -0.5 * params.pbcy_angstrom){
				dy += params.pbcy_angstrom;
			}
			if (dz > 0.5 * params.pbcz_angstrom){
				dz -= params.pbcz_angstrom;
			}
			if (dz < -0.5 * params.pbcz_angstrom){
				dz += params.pbcz_angstrom;
			}

			std::vector<double> force_comps = LJ_force(dx, dy, dz, params);
			system[i].fx += force_comps[0];
			system[i].fy += force_comps[1];
			system[i].fz += force_comps[2];

			system[j].fx -= force_comps[0];
			system[j].fy -= force_comps[1];
			system[j].fz -= force_comps[2];
		}
	}
}




void velocity_verlet(std::vector<Particle>& system, const Parameters& params, double dt){
	update_force_pbc(system, params);
	for (auto& p : system){
		p.vx += 0.5 * p.fx / params.mass_au * dt;
		p.vy += 0.5 * p.fy / params.mass_au * dt;
		p.vz += 0.5 * p.fz / params.mass_au * dt;

		p.x += p.vx * dt;
		p.y += p.vy * dt;
		p.z += p.vz * dt;
		
	}
	apply_pbc(system, params);
	update_force_pbc(system, params);
	for (auto& p : system){
		p.vx += 0.5 * p.fx / params.mass_au * dt;
		p.vy += 0.5 * p.fy / params.mass_au * dt;
		p.vz += 0.5 * p.fz / params.mass_au * dt;
	}
}	


int main(){
	Parameters params = read_config("config.txt");
	std::vector<Particle> system(params.N);
	initialize_system(system, params);
	std::ofstream fout(params.output);
	for (int i = 0; i < params.n_steps; ++i){
		velocity_verlet(system, params, params.dt);
		if (i % params.write_steps == 0){
			write_xyz(fout, system);
		}
	}
	fout.close();
	std::cout << "Finished simulation" << std::endl;
	return 0;
}



