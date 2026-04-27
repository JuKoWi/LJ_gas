#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <random>
#include "read_input.h"
#include "system.h"
#include "utils.h"
#include "energy_opt.h"


void System::initialize_random_positions(Parameters params){
	for (auto& p : particles){
		p.x = drand48() *  params.pbc_L_nm - 0.5 * params.pbc_L_nm;
		p.y = drand48() *  params.pbc_L_nm - 0.5 * params.pbc_L_nm;
		p.z = drand48() *  params.pbc_L_nm - 0.5 * params.pbc_L_nm;
	}
}

void System::initialize_positions_sc(Parameters params){
	std::cout << "Starting position initializion \n"; 
	int N_direction { static_cast<int>(std::ceil(std::cbrt(params.N))) };
	double dr { params.pbc_L_nm / N_direction };
	int count { 0 };
	for (int i=0; i<N_direction; ++i){
		for (int j=0; j<N_direction; ++j){
			for (int k=0; k<N_direction; ++k){
				if (count == params.N){
					return;
				}
				particles[count].x = i * dr;
				particles[count].y = j * dr;
				particles[count].z = k * dr;
				count += 1;
			}
		}
	}
}

void System::initialize_constant_velocity(){
	std::cout << "Initialize constant velocity \n";
	for (auto& p : particles){
		p.vx = 1;
	}
	return;
}

void System::initialize_rand_velocity(){
	std::cout << "Initialize random velocity \n";
	float vrange { 1e-3 };
	for (auto& p : particles){
		p.vx =drand48() *  vrange - 0.5 * vrange;
	}
	return;
}

void System::initialize_two_particles(Parameters params){
	particles[0].x = 0;
	particles[0].y = 0;
	particles[0].z = 0;

	particles[1].x = params.init_dist;
	particles[1].y = 0;
	particles[1].z = 0;
}

double System::get_temp(){
	return 2.0 * get_ekin() / (deg_of_freedom * PhysicalConstants::K_B);
}

double System::get_ekin(){
	double ekin = {0.0};
	for (const auto& p: particles){
		ekin += 0.5 * p.mass * (p.vx * p.vx + p.vy * p.vy + p.vz * p.vz);
	}
	return ekin;
}

// get center-of-mass momentum in g * nm/(mol * ps)
std::vector<double> System::get_cm_momentum(){
	double cm_momentum_x {0.0} ;
	double cm_momentum_y {0.0} ;
	double cm_momentum_z {0.0} ;
	for (const auto& p: particles){
		cm_momentum_x += p.mass * p.vx;
		cm_momentum_y += p.mass * p.vy;
		cm_momentum_z += p.mass * p.vz;
	}
	return {cm_momentum_x, cm_momentum_y, cm_momentum_z};
}

void System::initialize_temperature(double T_K){
	std::random_device rd;
	std::mt19937 gen(rd());

	for (auto& p: particles){
		double std_deviation = std::sqrt(PhysicalConstants::K_B * T_K / p.mass);
		std::normal_distribution dist(0.0, std_deviation);
		p.vx = dist(gen);
		p.vy = dist(gen);
		p.vz = dist(gen);
	}

}

void System::apply_vrescale(double T_K){
	if (T_K == 0){
		for (auto& p: particles){
			p.vx = 0.0;
			p.vy = 0.0;
			p.vz = 0.0;
		}
	}
	else{
		double scaling_factor { std::sqrt(T_K / get_temp()) };
		for (auto& p: particles){
			p.vx *= scaling_factor;
			p.vy *= scaling_factor;
			p.vz *= scaling_factor;
		}
	}
}

// double System::get_pressure(){
// 	double V {std::pow(box_size, 3)};
	
// }

// void System::do_simulation_step(){

// }

// Overall initialization to set up simulation
void System::initialize(const Parameters& params){
	box_size = params.pbc_L_nm;
	temperature = params.T_K;
	particles.resize(params.N);
	particles_at_verlet.resize(params.N);
	N = params.N;
	deg_of_freedom = 3 * N - 3 ;
	lj_sigma = params.sigma;
	square_lj_cutoff = params.sigma * params.sigma * std::pow(3,2);
	square_verlet_cutoff = params.sigma* params.sigma * std::pow(3 + 0.5, 2);
	


	for (auto& p: particles){
		p.mass = params.mass_au;
	}
	
	// initialize positions
	if (params.init_pos == "random"){
		initialize_random_positions(params);
	}
	else if (params.init_pos == "sc"){
		initialize_positions_sc(params);
	}
	else if (params.init_pos == "two"){
		if (params.N == 2){
			initialize_two_particles(params);
		}
		else{
			throw std::runtime_error("Wrong number of particles for this type of initial position");
		}
	}
	else {
		std::runtime_error("Undefined keyword for init_pos keyword");
	}

	// initialize velocities
	if (params.T_K > 0){
		initialize_temperature(params.T_K);
	}

	for (int i=0; i< N; i++){
		particles_at_verlet[i].x = particles[i].x;	
		particles_at_verlet[i].y = particles[i].y;	
		particles_at_verlet[i].z = particles[i].z;	
	}
}

bool System::verlet_update_required(){
	double square_verlet_skin {std::pow(std::sqrt(square_verlet_cutoff) - std::sqrt(square_lj_cutoff), 2)};
	for (int i=0; i<N; i++){
		if (particles[i].get_square_dist(particles_at_verlet[i], box_size) > square_verlet_skin/4){
			return true;
		}
	}
	return false;
}

void System::update_verlet(){
	if (verlet_update_required()){
		for (int i=0; i<N; i++){
			particles[i].verlet_list.clear();
			for (int j = 0; j < N; j++){
				if (i == j){
					continue;
				}
				if (particles[i].get_square_dist(particles[j], box_size) < square_verlet_cutoff){
					particles[i].verlet_list.push_back(j);
				}
			}
		}
		for (int i=0; i<N; i++){
		            particles_at_verlet[i].x = particles[i].x;
		            particles_at_verlet[i].y = particles[i].y;
		            particles_at_verlet[i].z = particles[i].z;
		}
	}
}

double Particle::get_square_dist(const Particle& particle, double box_size){
	double dx2 { std::pow(minimum_image_distance(x - particle.x, box_size), 2) };
	double dy2 { std::pow(minimum_image_distance(y - particle.y, box_size), 2) };
	double dz2 { std::pow(minimum_image_distance(z - particle.z, box_size), 2) };
	return dx2 + dy2 + dz2;
}
