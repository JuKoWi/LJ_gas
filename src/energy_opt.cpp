#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "energy_opt.h"
#include "system.h"
#include "utils.h"

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
			if (dx > 0.5 * params.pbc_L_angstrom){
				dx -= params.pbc_L_angstrom;
			}
			if (dx < -0.5 * params.pbc_L_angstrom){
				dx += params.pbc_L_angstrom;
			}
			if (dy > 0.5 * params.pbc_L_angstrom){
				dy -= params.pbc_L_angstrom;
			}
			if (dy < -0.5 * params.pbc_L_angstrom){
				dy += params.pbc_L_angstrom;
			}
			if (dz > 0.5 * params.pbc_L_angstrom){
				dz -= params.pbc_L_angstrom;
			}
			if (dz < -0.5 * params.pbc_L_angstrom){
				dz += params.pbc_L_angstrom;
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

void energy_opt(std::vector<Particle>& system, Parameters& params){
	double threshold = 0.0000001 * params.N;
	double total_force_sqr = 0;
	double eta = 0.001;
	int max_rep = 1000;
	std::ofstream fout_opt(params.output_opt);
	write_xyz(fout_opt, system);
	for (int i=0; i<max_rep; ++i){
		total_force_sqr = 0;
		update_force_pbc(system, params);
		for (auto& p: system){
			total_force_sqr += p.fx * p.fx;
			total_force_sqr += p.fy * p.fy;
			total_force_sqr += p.fz * p.fz;
		}
		std::cout << "Iteration: \n";
		std::cout << i << std::endl;
		std::cout << "Position particle 1 :\n";
		std::cout << system[0].x << std::endl;
		std::cout << "Total force: \n";
		std::cout << total_force_sqr << std::endl;
		if (total_force_sqr < threshold){
			std::cout << "Minimized energy below force-threshold \n";
			fout_opt.close();
			return;
		}
		for (auto& p:system){
			p.x += eta * p.fx;
			p.y += eta * p.fy;
			p.z += eta * p.fz;
		}
		update_pos_pbc(system, params);
		write_xyz(fout_opt, system);
	}
	fout_opt.close();
	std::cout << "Stopped energy minimization after 1000 iterations \n";
}

void update_pos_pbc(std::vector<Particle>& system, const Parameters& params){ // 0-centered simulation box
	for (auto& p : system){
		if (p.x > params.pbc_L_angstrom*0.5){
			p.x -= params.pbc_L_angstrom * std::floor(p.x / params.pbc_L_angstrom);
		}
		if (p.y > params.pbc_L_angstrom*0.5){
			p.y -= params.pbc_L_angstrom * std::floor(p.y / params.pbc_L_angstrom);
		}
		if (p.z > params.pbc_L_angstrom*0.5){
			p.z -= params.pbc_L_angstrom * std::floor(p.z / params.pbc_L_angstrom);
		}
		if (p.x < -params.pbc_L_angstrom*0.5){
			p.x -= params.pbc_L_angstrom * std::floor(p.x / params.pbc_L_angstrom);
		}
		if (p.y < -params.pbc_L_angstrom*0.5){
			p.y -= params.pbc_L_angstrom * std::floor(p.y / params.pbc_L_angstrom);
		}
		if (p.z < -params.pbc_L_angstrom*0.5){
			p.z -= params.pbc_L_angstrom * std::floor(p.z / params.pbc_L_angstrom);
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
	update_pos_pbc(system, params);
	update_force_pbc(system, params);
	for (auto& p : system){
		p.vx += 0.5 * p.fx / params.mass_au * dt;
		p.vy += 0.5 * p.fy / params.mass_au * dt;
		p.vz += 0.5 * p.fz / params.mass_au * dt;
	}
}	

double pair_potential_energy(double dr2, Parameters params){
	double U;
	U = 4 * params.epsilon * (std::pow(params.sigma*params.sigma/dr2, 6) - std::pow(params.sigma * params.sigma/dr2, 3));
	return U;
}

double potential_energy(const std::vector<Particle>& system, Parameters params){
	double dr2;
	double dx;
	double dy;
	double dz;
	double U = 0;
	for (int i=0; i<params.N; ++i){
		for (int j=i+1; j<params.N; ++j){
			dx = system[j].x - system[i].x;
			dy = system[j].y - system[i].y;
			dz = system[j].z - system[i].z;
			dr2 = dx * dx + dy * dy + dz * dz;
			U += pair_potential_energy(dr2, params);
		}
	}
	return U;
}



