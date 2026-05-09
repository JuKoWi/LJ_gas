#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "energy_opt.h"
#include "system.h"
#include "utils.h"

// Return components of LJ-force between two particles
// std::vector<double> LJ_force(double dx, double dy, double dz, const Parameters& params)
// {
//     double r2 {dx*dx + dy*dy + dz*dz};
// 	double r {std::sqrt(r2)};
//     if (r2 < 1e-12) return {0.0, 0.0, 0.0}; // TODO: use some maximal value, not 0

// 	double inv_r {1 / r};
// 	double f {24 * params.epsilon * (2 * std::pow(params.sigma * inv_r, 12) * inv_r - std::pow(params.sigma * inv_r, 6) * inv_r)}; //abs value force
// 	double dx_norm {dx * inv_r};
// 	double dy_norm {dy * inv_r};
// 	double dz_norm {dz * inv_r};
//     return {f * dx_norm, f * dy_norm, f * dz_norm};
// }

// std::vector<double> harmonic_force(double dx, double dy, double dz, const Parameters& params)
// {
//     double r2 { dx*dx + dy*dy + dz*dz };
// 	double r { std::sqrt(r2) };
// 	double inv_r { 1 / r };
//     if (r2 < 1e-12) return {0.0, 0.0, 0.0}; // TODO: use some maximal value, not 0

// 	double f {-2 * params.epsilon * (r - params.sigma)};
// 	double dx_norm { dx * inv_r };
// 	double dy_norm { dy * inv_r };
// 	double dz_norm { dz * inv_r };
//     return {f * dx_norm, f * dy_norm, f * dz_norm};
// }

// // distance components while considering PBC
// double minimum_image_distance(double dx, double L){
// 	if (dx > 0.5 * L){
// 		dx -= L;
// 	}
// 	if (dx < -0.5 * L){
// 		dx += L;
// 	}
// 	return dx;
// }

// // update forces taking pbc into account
// void update_force_pbc(System& system, const Parameters& params, bool verlet){
// 	for (auto& p : system.particles){
// 		p.fx = 0;
// 		p.fy = 0;
// 		p.fz = 0;
// 	}

// 	if (verlet){
// 		for (int i = 0; i<system.N; ++i){
// 			for (auto& j : system.particles[i].verlet_list){
// 				if (j<=i){
// 					continue;
// 				}
// 				double dx { minimum_image_distance(system.particles[j].x - system.particles[i].x, system.box_size) };
// 				double dy { minimum_image_distance(system.particles[j].y - system.particles[i].y, system.box_size) };
// 				double dz { minimum_image_distance(system.particles[j].z - system.particles[i].z, system.box_size) };

// 				if (dx*dx + dy*dy + dz*dz < system.square_lj_cutoff){
// 					std::vector<double> force_comps { LJ_force(dx, dy, dz, params) };
// 					// std::vector<double> force_comps = harmonic_force(dx, dy, dz, params);
// 					system.particles[i].fx -= force_comps[0];
// 					system.particles[i].fy -= force_comps[1];
// 					system.particles[i].fz -= force_comps[2];

// 					system.particles[j].fx += force_comps[0];
// 					system.particles[j].fy += force_comps[1];
// 					system.particles[j].fz += force_comps[2];
// 				}
// 			}
// 		}
// 	}
// 	else{
// 		for (int i = 0; i<system.N; ++i){
// 			for (int j = i + 1; j < system.N; ++j){ // only iterate over pairs
// 				double dx { minimum_image_distance(system.particles[j].x - system.particles[i].x, system.box_size) };
// 				double dy { minimum_image_distance(system.particles[j].y - system.particles[i].y, system.box_size) };
// 				double dz { minimum_image_distance(system.particles[j].z - system.particles[i].z, system.box_size) };

// 				if (dx*dx + dy*dy + dz*dz < system.square_lj_cutoff){
// 					std::vector<double> force_comps { LJ_force(dx, dy, dz, params) };
// 					// std::vector<double> force_comps = harmonic_force(dx, dy, dz, params);
// 					system.particles[i].fx -= force_comps[0];
// 					system.particles[i].fy -= force_comps[1];
// 					system.particles[i].fz -= force_comps[2];

// 					system.particles[j].fx += force_comps[0];
// 					system.particles[j].fy += force_comps[1];
// 					system.particles[j].fz += force_comps[2];
// 				}

// 			}
// 		}
// 	}
// }

// // steepest descend energy optimization
// void energy_opt(System& system, const Parameters& params){
// 	double threshold { 1e-9 * params.epsilon * params.N };
// 	double total_force_sqr { 0 };
// 	double eta { 1 };
// 	double maxstep { 0.01 }; //nm
// 	int max_rep { 1000 };
// 	std::ofstream fout_opt(params.output_opt);
// 	write_xyz_frame(fout_opt, system);
// 	double Epot_before {};
// 	double Epot_after {};

// 	for (int i=0; i<max_rep; ++i){
// 		total_force_sqr = 0;
// 		update_force_pbc(system, params, false);
// 		for (auto& p: system.particles){
// 			total_force_sqr += p.fx * p.fx;
// 			total_force_sqr += p.fy * p.fy;
// 			total_force_sqr += p.fz * p.fz;
// 		}
// 		if (i % 100 == 0){
// 			std::cout << "Iteration: " << i << " Total force: " << total_force_sqr << '\n';
// 		}
// 		if (total_force_sqr < threshold){
// 			std::cout << "Minimized energy below force-threshold \n";
// 			fout_opt.close();
// 			return;
// 		}

// 		Epot_before = potential_energy(system, params);
// 		for (auto& p: system.particles){
// 			double dx { eta * p.fx };
// 			p.x += (std::abs(dx) < maxstep) ? dx : std::copysign(maxstep, dx);
// 			double dy {eta * p.fy };
// 			p.y += (std::abs(dy) < maxstep) ? dy : std::copysign(maxstep, dy);
// 			double dz { eta * p.fz };
// 			p.z += (std::abs(dz) < maxstep) ? dz : std::copysign(maxstep, dz);
// 		}
// 		update_pos_pbc(system);
// 		write_xyz_frame(fout_opt, system);
// 		Epot_after = potential_energy(system, params);
// 		if (Epot_after > Epot_before){
// 			eta *= 0.5;
// 		}
// 		else{
// 			eta *= 1.1;
// 		}
// 	}
// 	fout_opt.close();
// 	std::cout << "Stopped energy minimization after 1000 iterations \n";
// }

// // map positions outside simulation box back to simulation box
// void update_pos_pbc(System& system){ // 0-centered simulation box
// 	double L { system.box_size};
// 	double half_box { L * 0.5 };

// 	for (auto& p : system.particles){
// 		p.x -= L * std::floor((p.x + half_box) / L);
// 		p.y -= L * std::floor((p.y + half_box) / L);
// 		p.z -= L * std::floor((p.z + half_box) / L);
// 	}
// }

// // do a velocity verlet step
// void velocity_verlet(System& system, const Parameters& params, double dt){
// 	system.update_verlet();
// 	update_force_pbc(system, params, true);
// 	for (auto& p : system.particles){
// 		p.vx += 0.5 * p.fx / params.mass_au * dt;
// 		p.vy += 0.5 * p.fy / params.mass_au * dt;
// 		p.vz += 0.5 * p.fz / params.mass_au * dt;

// 		p.x += p.vx * dt;
// 		p.y += p.vy * dt;
// 		p.z += p.vz * dt;
// 	}
// 	update_pos_pbc(system);
// 	update_force_pbc(system, params, true);
// 	for (auto& p : system.particles){
// 		p.vx += 0.5 * p.fx / params.mass_au * dt;
// 		p.vy += 0.5 * p.fy / params.mass_au * dt;
// 		p.vz += 0.5 * p.fz / params.mass_au * dt;
// 	}
// }	

// // LJ potential contribution for one pair in kj/mol
// double pair_potential_energy(double dr2, const Parameters params){
// 	double sqr_sigma { params.sigma * params.sigma };
// 	return params.epsilon + 4 * params.epsilon * (std::pow(sqr_sigma / dr2, 6) - std::pow(sqr_sigma / dr2, 3));
// }

// double pair_potential_harmonic(double dr2, Parameters params){
// 	double r { std::sqrt(dr2) };
// 	return params.epsilon * std::pow(r - params.sigma, 2);
// }

// // total potential energy in kj/mol
// double potential_energy(const System& system, const Parameters params){
// 	double dr2 {};
// 	double dx {};
// 	double dy {};
// 	double dz {};
// 	double U { 0 };
// 	for (int i=0; i<system.N; ++i){
// 		for (int j=i+1; j<system.N; ++j){
// 			dx = minimum_image_distance(system.particles[j].x - system.particles[i].x, system.box_size);
// 			dy = minimum_image_distance(system.particles[j].y - system.particles[i].y, system.box_size);
// 			dz = minimum_image_distance(system.particles[j].z - system.particles[i].z, system.box_size);
// 			dr2 = dx * dx + dy * dy + dz * dz;
// 			U += pair_potential_energy(dr2, params);
// 			// U += pair_potential_harmonic(dr2, params);
// 		}
// 	}
// 	return U;
// }

// // total kinetic energy in kj/mol
// double kinetic_energy(const System& system){
// 	double Ekin { 0 };
// 	for (auto& p : system.particles){
// 		Ekin += p.mass/2 * p.vx * p.vx; // component wise contribution in kj/mol
// 		Ekin += p.mass/2 * p.vy * p.vy; // component wise contribution in kj/mol
// 		Ekin += p.mass/2 * p.vz * p.vz; // component wise contribution in kj/mol
// 	}
// 	return Ekin;
// }

// // find the maximal force on all particles

