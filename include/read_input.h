#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <vector>
#include <unordered_map>
#include "system.h"

enum class Thermostat{
	None,
	VRescale,
	Berendsen,
	NoseHoover,
};

// struct Parameters{
// 	int N;
// 	int n_steps;
// 	double dt;
// 	double epsilon;
// 	double sigma;
// 	double mass_au;
// 	double pbc_L_nm;
// 	double T_K;
// 	int optimization;
// 	int write_steps;
// 	std::string output;
// 	std::string output_opt;
// 	std::string init_pos;
// 	double init_dist; 
// 	std::string ensemble_type;
// 	Thermostat thermostat;
// };

struct JobParameters{
	int n_steps {};
	double dt {};
	double T_K {};
	int optimization {};
	int write_steps {};
	std::string outfile_traj_lammps {};
	std::string outfile_opt_xyz {};
	std::string outfile_energy {};
	std::string infile_geom_json {};
	std::string ff_dir {};
	std::string ensemble_type {};
	Thermostat thermostat {};
};

JobParameters parse_job(const std::string& filename);
Thermostat parse_thermostat(const std::string& key_thermostat);
std::vector<std::string> get_forcefield_files(const std::string& ff_dir);

struct AtomType; // why is this necessary to declare here?
struct BondType; // why is this necessary to declare here?
std::unordered_map<std::string, AtomType> parse_atom_types(const std::vector<std::string>& file_list);

std::string make_bond_key(std::string string_a, std::string string_b);
std::vector<BondType> parse_bond_types(const std::vector<std::string>& file_list, const std::unordered_map<std::string, AtomType>& atom_types);
bool bondtype_exists(std::vector<BondType> bond_types, std::string name_of_type);
int get_type_idx(std::vector<BondType> bond_types, std::string name_of_type);

#endif
