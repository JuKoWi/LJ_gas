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
struct AngleType; // why is this necessary to declare here?
std::unordered_map<std::string, AtomType> parse_atom_types(const std::vector<std::string>& file_list);

std::string make_bond_key(std::string string_a, std::string string_b);
std::string make_angle_key(std::string string_a, std::string string_b, std::string string_c);
std::vector<BondType> parse_bond_types(const std::vector<std::string>& file_list, const std::unordered_map<std::string, AtomType>& atom_types);
std::vector<AngleType> parse_angle_types(const std::vector<std::string>& file_list, const std::unordered_map<std::string, AtomType>& atom_types);

template <typename T>
bool type_exists(std::vector<T> types, std::string name_of_type){
	for (auto t: types){
		if (t.name == name_of_type){
			return true;
		}
	}
	return false;
}

template<typename T>
int get_type_idx(const std::vector<T>& types, const std::string& name_of_type){
	for (size_t i=0; i < types.size(); ++i){
		if (types[i].name == name_of_type){
			return static_cast<int>(i);
		}
	}
	return -1;
}

#endif
