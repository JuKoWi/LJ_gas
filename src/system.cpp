#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <random>
#include <cassert>
#include "json.hpp"
#include "read_input.h"
#include "system.h"
#include "utils.h"
#include "energy_opt.h"

using json = nlohmann::json;


void System::initialize_atoms(const std::string& geom_file, std::unordered_map<std::string, AtomType> type_dict){
    std::cout << "Initialize atom positions and connectivity \n";
    atom_type_dict = type_dict;
	std::ifstream file(geom_file);
	if (!file){
		throw std::runtime_error("Cannot open geometry file!");
	}
	json data;
	file >> data;

    for (auto L: data["box"]){
        box_size.push_back(L);
    }

    for (auto atom: data["atoms"]){
        atom_types.push_back(atom["type"]);
        x.push_back(atom["pos"][0]);
        y.push_back(atom["pos"][1]);
        z.push_back(atom["pos"][2]);
        vx.push_back(0.0);
        vy.push_back(0.0);
        vz.push_back(0.0);
        fx.push_back(0.0);
        fy.push_back(0.0);
        fz.push_back(0.0);
        sigma_nm.push_back(atom_type_dict.at(atom["type"]).sigma);
        charges.push_back(atom_type_dict.at(atom["type"]).charge);
        epsilon_kjmol.push_back(atom_type_dict.at(atom["type"]).epsilon);
        masses_gmol.push_back(atom_type_dict.at(atom["type"]).mass);
    }
    N = atom_types.size();
    deg_of_freedom = 3 * N - 3;
    set_lj_pairs();

    assert(atom_types.size() == x.size());
    assert(x.size() == y.size());
    assert(y.size() == z.size());
    assert(x.size() == fx.size());
    assert(x.size() == fy.size());
    assert(x.size() == fz.size());
    assert(pair_epsilon_kjmol.size() == x.size()*x.size());
    assert(pair_sigma_nm.size() == x.size()*x.size());
}

void System::set_lj_pairs(){
    pair_epsilon_kjmol.resize(N * N);
    pair_sigma_nm.resize(N * N);
    for (int i=0; i<N; ++i){
        for (int j=0; j<N; ++j){
            pair_sigma_nm[i*N + j] = 0.5 * (sigma_nm[i] + sigma_nm[j]);
            pair_epsilon_kjmol[i*N + j] = std::sqrt(epsilon_kjmol[i] * epsilon_kjmol[j]);
        }
    }
}

void System::do_verlet_step(double dt){
    //TODO: update verlet
    update_forces();
    for (int i=0; i<N; i++){
        vx[i] = 0.5 * fx[i] / masses_gmol[i] * dt;
        vy[i] = 0.5 * fy[i] / masses_gmol[i] * dt;
        vz[i] = 0.5 * fz[i] / masses_gmol[i] * dt;

        x[i] += vx[i] * dt;
        y[i] += vy[i] * dt;
        z[i] += vz[i] * dt;
    }
    update_positions();
    update_forces();
    for (int i=0; i<N; i++){
        vx[i] = 0.5 * fx[i] / masses_gmol[i] * dt;
        vy[i] = 0.5 * fy[i] / masses_gmol[i] * dt;
        vz[i] = 0.5 * fz[i] / masses_gmol[i] * dt;
    }
}

void System::do_simulation_step(double dt, int step_count){
    do_verlet_step(dt);
}

void System::update_forces(){
    std::fill(fx.begin(), fx.end(), 0.0);
    std::fill(fy.begin(), fy.end(), 0.0);
    std::fill(fz.begin(), fz.end(), 0.0);

    update_lj_forces();
}

void System::update_lj_forces(){
    for (int i=0; i<N; ++i){
        for (int j=i+1; j<N; ++j){
            double dx {minimum_image_distance(x[j] - x[i], box_size[0])};
            double dy {minimum_image_distance(y[j] - y[i], box_size[1])};
            double dz {minimum_image_distance(z[j] - z[i], box_size[2])};
            if (dx*dx + dy*dy + dz*dz < std::pow(square_lj_cutoff * pair_sigma_nm[i*N + j], 2)){
                Vec3 f_lj {LJ_force(dx, dy, dz, pair_sigma_nm[i*N + j], pair_epsilon_kjmol[i*N +j])};
				fx[i] -= f_lj.x;
				fy[i] -= f_lj.y;
				fz[i] -= f_lj.z;

				fx[j] += f_lj.x;
				fy[j] += f_lj.y;
				fz[j] += f_lj.z;
            }
        }
    }
}

// // distance components while considering PBC
double minimum_image_distance(double dx, double L){
	if (dx > 0.5 * L){
		dx -= L;
	}
	if (dx < -0.5 * L){
		dx += L;
	}
	return dx;
}

// Return components of LJ-force between two particles
inline Vec3 LJ_force(double dx, double dy, double dz, double sigma, double epsilon) {
    double r2 {dx*dx + dy*dy + dz*dz};
    Vec3 f_vec {};
    if (r2 < 1e-12) return f_vec; // TODO: use some maximal value, not 0
	double r {std::sqrt(r2)};

	double inv_r {1 / r};
	double f {24 * epsilon * (2 * std::pow(sigma * inv_r, 12) * inv_r - std::pow(sigma * inv_r, 6) * inv_r)}; //abs value force
	double dx_norm {dx * inv_r};
	double dy_norm {dy * inv_r};
	double dz_norm {dz * inv_r};
    f_vec.x = f * dx_norm;
    f_vec.y = f * dy_norm; 
    f_vec.z = f * dz_norm;
    return f_vec;
}

// map positions outside simulation box back to simulation box
void System::update_positions(){ // 0-centered simulation box
    double half_box_x {box_size[0]/2};
    double half_box_y {box_size[1]/2};
    double half_box_z {box_size[2]/2};
	for (int i=0; i<N; ++i){
		x[i] -= box_size[0] * std::floor((x[i] + half_box_x)/ box_size[0]);
		y[i] -= box_size[1] * std::floor((y[i] + half_box_y)/ box_size[1]);
		z[i] -= box_size[2] * std::floor((z[i] + half_box_z)/ box_size[2]);
	}
}

void System::run_simulation(JobParameters params){
	std::ofstream fout_traj(params.outfile_traj_lammps);
	double time {};
	for (int i = 0; i < params.n_steps; ++i){
		do_simulation_step(params.dt, i);
		if (i % params.write_steps == 0){
			write_lammpsdump_frame(fout_traj, i);
			// write_xyz_frame(fout_traj, system);
			// time = i * params.dt ;// time in ps
			// Ekin = kinetic_energy(system);
			// Epot = potential_energy(system, params);
			// write_frame_diagnostics(fout_energy, time, Ekin, Epot, system.get_temp());
		}
	}
}


// void System::initialize_random_positions(Parameters params){
// 	for (auto& p : particles){
// 		p.x = drand48() *  params.pbc_L_nm - 0.5 * params.pbc_L_nm;
// 		p.y = drand48() *  params.pbc_L_nm - 0.5 * params.pbc_L_nm;
// 		p.z = drand48() *  params.pbc_L_nm - 0.5 * params.pbc_L_nm;
// 	}
// }

// void System::initialize_positions_sc(Parameters params){
// 	std::cout << "Starting position initializion \n"; 
// 	int N_direction { static_cast<int>(std::ceil(std::cbrt(params.N))) };
// 	double dr { params.pbc_L_nm / N_direction };
// 	int count { 0 };
// 	for (int i=0; i<N_direction; ++i){
// 		for (int j=0; j<N_direction; ++j){
// 			for (int k=0; k<N_direction; ++k){
// 				if (count == params.N){
// 					return;
// 				}
// 				particles[count].x = i * dr;
// 				particles[count].y = j * dr;
// 				particles[count].z = k * dr;
// 				count += 1;
// 			}
// 		}
// 	}
// }

// void System::initialize_constant_velocity(){
// 	std::cout << "Initialize constant velocity \n";
// 	for (auto& p : particles){
// 		p.vx = 1;
// 	}
// 	return;
// }

// void System::initialize_rand_velocity(){
// 	std::cout << "Initialize random velocity \n";
// 	float vrange { 1e-3 };
// 	for (auto& p : particles){
// 		p.vx =drand48() *  vrange - 0.5 * vrange;
// 	}
// 	return;
// }

// void System::initialize_two_particles(Parameters params){
// 	particles[0].x = 0;
// 	particles[0].y = 0;
// 	particles[0].z = 0;

// 	particles[1].x = params.init_dist;
// 	particles[1].y = 0;
// 	particles[1].z = 0;
// }

// double System::get_temp(){
// 	return 2.0 * get_ekin() / (deg_of_freedom * PhysicalConstants::K_B);
// }

// double System::get_ekin(){
// 	double ekin = {0.0};
// 	for (const auto& p: particles){
// 		ekin += 0.5 * p.mass * (p.vx * p.vx + p.vy * p.vy + p.vz * p.vz);
// 	}
// 	return ekin;
// }

// // get center-of-mass momentum in g * nm/(mol * ps)
// std::vector<double> System::get_cm_momentum(){
// 	double cm_momentum_x {0.0} ;
// 	double cm_momentum_y {0.0} ;
// 	double cm_momentum_z {0.0} ;
// 	for (const auto& p: particles){
// 		cm_momentum_x += p.mass * p.vx;
// 		cm_momentum_y += p.mass * p.vy;
// 		cm_momentum_z += p.mass * p.vz;
// 	}
// 	return {cm_momentum_x, cm_momentum_y, cm_momentum_z};
// }

// void System::initialize_temperature(double T_K){
// 	std::random_device rd;
// 	std::mt19937 gen(rd());

// 	for (auto& p: particles){
// 		double std_deviation = std::sqrt(PhysicalConstants::K_B * T_K / p.mass);
// 		std::normal_distribution dist(0.0, std_deviation);
// 		p.vx = dist(gen);
// 		p.vy = dist(gen);
// 		p.vz = dist(gen);
// 	}

// }

// void System::apply_vrescale(double T_K){
// 	if (T_K == 0){
// 		for (auto& p: particles){
// 			p.vx = 0.0;
// 			p.vy = 0.0;
// 			p.vz = 0.0;
// 		}
// 	}
// 	else{
// 		double scaling_factor { std::sqrt(T_K / get_temp()) };
// 		for (auto& p: particles){
// 			p.vx *= scaling_factor;
// 			p.vy *= scaling_factor;
// 			p.vz *= scaling_factor;
// 		}
// 	}
// }

// // double System::get_pressure(){
// // 	double V {std::pow(box_size, 3)};
	
// // }

// // void System::do_simulation_step(){

// // }

// void System::initialize(const Parameters& params){
// 	box_size = params.pbc_L_nm;
// 	temperature = params.T_K;

// 	x.resize(params.N);
// 	y.resize(params.N);
// 	z.resize(params.N);

// 	x_at_verlet.resize(params.N);
// 	y_at_verlet.resize(params.N);
// 	z_at_verlet.resize(params.N);

// 	vx.resize(params.N);
// 	vy.resize(params.N);
// 	vz.resize(params.N);

// 	fx.resize(params.N);
// 	fy.resize(params.N);
// 	fz.resize(params.N);

// 	masses_gmol.resize(params.N);
// 	charges.resize(params.N);

// 	N = params.N;
// 	deg_of_freedom = 3 * N - 3 ;


// }

// // Overall initialization to set up simulation
// void System::initialize(const Parameters& params){
// 	box_size = params.pbc_L_nm;
// 	temperature = params.T_K;
// 	particles.resize(params.N);
// 	particles_at_verlet.resize(params.N);
// 	N = params.N;
// 	deg_of_freedom = 3 * N - 3 ;
// 	lj_sigma = params.sigma;
// 	square_lj_cutoff = params.sigma * params.sigma * std::pow(3,2);
// 	square_verlet_cutoff = params.sigma* params.sigma * std::pow(3 + 0.5, 2);
	


// 	for (auto& p: particles){
// 		p.mass = params.mass_au;
// 	}
	
// 	// initialize positions
// 	if (params.init_pos == "random"){
// 		initialize_random_positions(params);
// 	}
// 	else if (params.init_pos == "sc"){
// 		initialize_positions_sc(params);
// 	}
// 	else if (params.init_pos == "two"){
// 		if (params.N == 2){
// 			initialize_two_particles(params);
// 		}
// 		else{
// 			throw std::runtime_error("Wrong number of particles for this type of initial position");
// 		}
// 	}
// 	else {
// 		std::runtime_error("Undefined keyword for init_pos keyword");
// 	}

// 	// initialize velocities
// 	if (params.T_K > 0){
// 		initialize_temperature(params.T_K);
// 	}

// 	for (int i=0; i< N; i++){
// 		particles_at_verlet[i].x = particles[i].x;	
// 		particles_at_verlet[i].y = particles[i].y;	
// 		particles_at_verlet[i].z = particles[i].z;	
// 	}
// }

// bool System::verlet_update_required(){
// 	double square_verlet_skin {std::pow(std::sqrt(square_verlet_cutoff) - std::sqrt(square_lj_cutoff), 2)};
// 	for (int i=0; i<N; i++){
// 		if (particles[i].get_square_dist(particles_at_verlet[i], box_size) > square_verlet_skin/4){
// 			return true;
// 		}
// 	}
// 	return false;
// }

// void System::update_verlet(){
// 	if (verlet_update_required()){
// 		for (int i=0; i<N; i++){
// 			particles[i].verlet_list.clear();
// 			for (int j = 0; j < N; j++){
// 				if (i == j){
// 					continue;
// 				}
// 				if (particles[i].get_square_dist(particles[j], box_size) < square_verlet_cutoff){
// 					particles[i].verlet_list.push_back(j);
// 				}
// 			}
// 		}
// 		for (int i=0; i<N; i++){
// 		            particles_at_verlet[i].x = particles[i].x;
// 		            particles_at_verlet[i].y = particles[i].y;
// 		            particles_at_verlet[i].z = particles[i].z;
// 		}
// 	}
// }

