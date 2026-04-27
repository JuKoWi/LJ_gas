#ifndef SYSTEM_H
#define SYSTEM_H

#include "read_input.h"
#include <vector>

struct System;

// hold particle-wise information
struct Particle{
	double x {}; // nm
	double y {};
	double z {};  

	double vx {}; // nm/ps
	double vy {};
	double vz {};  

	double fx {}; // kj/(mol nm)
	double fy {};
	double fz {};  

	double mass {};
	std::vector<int> verlet_list {};

	double get_square_dist(const Particle& particle, double box_size);
	
};


// hold all instantaneous information (particles, box size, temp,...)
struct System{
	std::vector<Particle> particles {};
	std::vector<Particle> particles_at_verlet {}; //keep track of positions at the time of Verlet list update
	double box_size {};
	double temperature {};
	int N {};
	int deg_of_freedom {};
	double lj_sigma {};
	double square_lj_cutoff {};
	double square_verlet_cutoff {};


	// Overall initialization
	void initialize(const Parameters& params);

	// Steps during initialization
	void initialize_random_positions(Parameters params);
	void initialize_positions_sc(Parameters params);
	void initialize_constant_velocity();
	void initialize_rand_velocity();
	void initialize_positions_file();
	void initialize_two_particles(Parameters params);
	void initialize_temperature(double T_K);

	//write configuration to file

	//get thermodynamics
	double get_ekin();
	double get_temp();
	std::vector<double> get_cm_momentum();
	double get_pressure();

	//thermostats
	void apply_vrescale(double T_K);

	
	void do_verlet_step();

	// do all actions required in one simulation step
	void do_simulation_step();
	bool verlet_update_required();
	void update_verlet();


};

#endif
