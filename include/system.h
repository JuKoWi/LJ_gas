#ifndef SYSTEM_H
#define SYSTEM_H

#include "read_input.h"
#include <vector>

struct System;
struct JobParameters;

struct Bonds{
	std::vector<int> i {};
	std::vector<int> j {};
	std::vector<int> type {};
};

struct BondType{
	std::string name {};
	double r0 {};
	double k {};

	void print(){
		std::cout << "BONDTYPE: " << name << '\n';
		std::cout << "\tr0: " << r0 << '\n';
		std::cout << "\tk: " << k << '\n';
	}
};

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
	double temperature_goal {};
	int N {};
	int deg_of_freedom {};
	double lj_cutoff_nm {};
	double verlet_skin_nm {0.2};
	double verlet_cutoff_nm {};

	std::vector<double> x {}; // nm
	std::vector<double> y {}; // nm
	std::vector<double> z {}; // nm

	std::vector<double> x_at_verlet {}; //keep track of positions at the time of Verlet list update
	std::vector<double> y_at_verlet {};
	std::vector<double> z_at_verlet {};
	std::vector<std::vector<int>> verlet_list {}; // neighbor list, implement flat list later

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

	std::vector<BondType> bond_types {};
	Bonds bonds {};

	// Steps during initialization
	void initialize(JobParameters job_params, std::unordered_map<std::string, AtomType> atom_types, std::vector<BondType> bond_types);
	void initialize_atoms(const std::string& geom_file, std::unordered_map<std::string, AtomType> type_dict);
	void initialize_bonds(const std::string& geom_file, std::vector<BondType> bond_types);
	void set_lj_pairs();
	void initialize_temperature(double T_K);
	void energy_opt(std::string opt_out);

	//main loop
	void do_verlet_step(double dt);
	void update_forces();
	void update_lj_forces(bool verlet);
	void update_bond_forces();
	void update_positions();

	// do all actions required in one simulation step
	void do_simulation_step(double dt);
	void run_simulation(JobParameters params);
	bool verlet_update_required();
	void update_verlet();
	double get_square_distance(int i, int j);

	//get thermodynamics
	double get_temp();
	Vec3 get_cm_momentum();
	// double get_pressure();
	double get_kinetic_energy();
	double lj_pair_potential(double dr_sqare, int id1, int id2);
	double bond_pair_potential(double dr, double r0, double k);
	double get_total_lj_potential();
	double get_total_bond_potential();
	double get_total_potential();
	double get_total_energy();

	//thermostats
	void apply_vrescale(double T_K);


	//write configuration to file
	void write_xyz_frame(std::ofstream& fout);
	void write_lammpsdump_frame(std::ofstream& fout, int step_count);
	void write_frame_diagnostics(std::ofstream& fout, double time, double Ekin, double Epot, double T_K);

};

double minimum_image_distance(double dx, double L);
inline Vec3 lj_force(double dx, double dy, double dz, double sigma, double epsilon);
inline Vec3 bond_force(double dx, double dy, double dz, double r0, double k);

#endif
