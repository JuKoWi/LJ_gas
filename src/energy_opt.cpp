#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "energy_opt.h"
#include "system.h"
#include "utils.h"

// Return LJ-force between two particles
std::vector<double> LJ_force(double dx, double dy, double dz, const Parameters& params)
{
    double r2 = dx*dx + dy*dy + dz*dz;
	double r = std::sqrt(r2);
    if (r2 < 1e-12) return {0.0, 0.0, 0.0}; // TODO: use some maximal value, not 0

	double inv_r = 1 / r;
	double f = 24 * params.epsilon * (2 * std::pow(params.sigma * inv_r, 12) * inv_r - std::pow(params.sigma * inv_r, 6) * inv_r);
    return {f * dx, f * dy, f * dz};
}

// update forces taking pbc into account
void update_force_pbc(std::vector<Particle>& system, const Parameters& params){
	double square_cutoff = 9; // interaction cutoff for the square distance 
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
			if (dx > 0.5 * params.pbc_L_nm){
				dx -= params.pbc_L_nm;
			}
			if (dx < -0.5 * params.pbc_L_nm){
				dx += params.pbc_L_nm;
			}
			if (dy > 0.5 * params.pbc_L_nm){
				dy -= params.pbc_L_nm;
			}
			if (dy < -0.5 * params.pbc_L_nm){
				dy += params.pbc_L_nm;
			}
			if (dz > 0.5 * params.pbc_L_nm){
				dz -= params.pbc_L_nm;
			}
			if (dz < -0.5 * params.pbc_L_nm){
				dz += params.pbc_L_nm;
			}

			if (dx*dx + dy*dy + dz*dz < square_cutoff){
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
}

void energy_opt(std::vector<Particle>& system, Parameters& params){
	double threshold = 1e-3 * params.N;
	double total_force_sqr = 0;
	double eta = 1e-3;
	double maxstep = 0.01; //nm
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
		if (i % 100 == 0){
			std::cout << "Iteration: " << i << " Total force: " << total_force_sqr << std::endl;
		}
		if (total_force_sqr < threshold){
			std::cout << "Minimized energy below force-threshold \n";
			fout_opt.close();
			return;
		}
		double Epot_before = potential_energy(system, params);
		for (auto& p:system){
			double dx = eta * p.fx;
			p.x += (std::abs(dx) < maxstep) ? dx : std::copysign(maxstep, dx);
			double dy = eta * p.fy;
			p.y += (std::abs(dy) < maxstep) ? dy : std::copysign(maxstep, dy);
			double dz = eta * p.fz;
			p.z += (std::abs(dz) < maxstep) ? dz : std::copysign(maxstep, dz);
		}
		update_pos_pbc(system, params);
		write_xyz(fout_opt, system);
		double Epot_after = potential_energy(system, params);
		if (Epot_after > Epot_before){
			eta *= 0.5;
		}
		else{
			eta *= 1.05;
		}
	}
	fout_opt.close();
	std::cout << "Stopped energy minimization after 1000 iterations \n";
}

// map positions outside simulation box back to simulation box
void update_pos_pbc(std::vector<Particle>& system, const Parameters& params){ // 0-centered simulation box
	double L = params.pbc_L_nm;
	double half_box = L * 0.5;
	for (auto& p : system){
		p.x -= L * std::floor((p.x + half_box) / L);
		p.y -= L * std::floor((p.y + half_box) / L);
		p.z -= L * std::floor((p.z + half_box) / L);
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

// LJ potential contribution for one pair in kj/mol
double pair_potential_energy(double dr2, Parameters params){
	double sqr_sigma = params.sigma * params.sigma;
	double U = 4 * params.epsilon * (std::pow(sqr_sigma / dr2, 6) - std::pow(sqr_sigma / dr2, 3));
	return U;
}

// total potential energy in kj/mol
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

// total kinetic energy in kj/mol
double kinetic_energy(const std::vector<Particle>& system, Parameters params){
	double Ekin = 0;
	for (auto& p : system){
		Ekin += params.mass_au/2 * p.vx * p.vx; // component wise contribution in kj/mol
		Ekin += params.mass_au/2 * p.vy * p.vy; // component wise contribution in kj/mol
		Ekin += params.mass_au/2 * p.vz * p.vz; // component wise contribution in kj/mol
	}
	return Ekin;
}

// find the maximal force on all particles

