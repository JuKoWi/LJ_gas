#ifndef SYSTEM_H
#define SYSTEM_H

#include "read_input.h"
#include <vector>

struct System;
struct JobParameters;

struct Vec3{
	double x {};
	double y {};
	double z {};
};

struct AtomType{
	std::string name {};
	std::string element {};
	double mass {};
	double charge {};
	double sigma {};
	double epsilon {};	

	void print(){
		std::cout << "TYPE: " << name << '\n';
		std::cout << "\tMass: " << mass << '\n';
		std::cout << "\tSigma: " << sigma << '\n';
		std::cout << "\tEpsilon: " << epsilon << '\n';
	}
};


// hold all instantaneous information (particles, box size, temp,...)
struct System{
	std::vector<double> box_size {};
	double temperature {};
	int N {};
	int deg_of_freedom {};
	double square_lj_cutoff {2.5};
	double square_verlet_cutoff {3.5};

	std::vector<double> x {}; // nm
	std::vector<double> y {}; // nm
	std::vector<double> z {}; // nm

	std::vector<double> x_at_verlet {}; //keep track of positions at the time of Verlet list update
	std::vector<double> y_at_verlet {};
	std::vector<double> z_at_verlet {};

	std::vector<double> vx {}; // nm/ps
	std::vector<double> vy {}; // nm/ps
	std::vector<double> vz {}; // nm/ps

	std::vector<double> fx {}; // kj/(mol nm)
	std::vector<double> fy {}; // kj/(mol nm)
	std::vector<double> fz {}; // kj/(mol nm)

	std::vector<double> masses_gmol {};
	std::vector<double> charges {};
	std::vector<double> epsilon_kjmol {};
	std::vector<double> sigma_nm {};
	std::vector<double> pair_sigma_nm {};
	std::vector<double> pair_epsilon_kjmol {};
	std::unordered_map<std::string, AtomType> atom_type_dict {};
	std::vector<std::string> atom_types {};

	std::vector<std::vector<int>> verlet_list {}; // neighbor list, implement flat list later

	double Ekin {};
	double Epot {};

	// Steps during initialization
	void initialize_atoms(const std::string& geom_file, std::unordered_map<std::string, AtomType> type_list);
	void set_lj_pairs();
	// void initialize_positions_sc(Parameters params);
	// void initialize_constant_velocity();
	// void initialize_rand_velocity();
	// void initialize_positions_file();
	// void initialize_two_particles(Parameters params);
	// void initialize_temperature(double T_K);

	// Overall initialization
	// void initialize(const Parameters& params);

	//write configuration to file

	//get thermodynamics
	// double get_ekin();
	// double get_temp();
	// std::vector<double> get_cm_momentum();
	// double get_pressure();

	//thermostats
	// void apply_vrescale(double T_K);

	void do_verlet_step(double dt);
	void update_forces();
	void update_lj_forces();
	void update_positions();

	// do all actions required in one simulation step
	void do_simulation_step(double dt, int step_count);
	void run_simulation(JobParameters params);
	// bool verlet_update_required();
	// void update_verlet();

	void write_xyz_frame(std::ofstream& fout);
	void write_lammpsdump_frame(std::ofstream& fout, int step_count);
};

double minimum_image_distance(double dx, double L);
inline Vec3 LJ_force(double dx, double dy, double dz, double sigma, double epsilon);

#endif
