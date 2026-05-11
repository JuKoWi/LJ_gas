#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <filesystem>
#include "read_input.h"
#include "json.hpp"
#include "system.h"

using json = nlohmann::json;

JobParameters parse_job(const std::string& filename){
	std::cout << "Start reading job file \n";
	std::ifstream file(filename);
	if (!file){
		throw std::runtime_error("Cannot open config file!");
	}

	json data;
	file >> data;

	JobParameters job {};
	job.dt = data["dt"];
	job.n_steps = data["n_steps"];
	job.T_K = data["T_K"];
	job.write_steps = data["write_interval"];
	job.outfile_traj_lammps = data["paths"]["out_traj"].get<std::string>();
	job.outfile_opt_xyz = data["paths"]["out_opt"].get<std::string>();
	job.outfile_energy = data["paths"]["out_energy"].get<std::string>();
	job.infile_geom_json = data["paths"]["geom_json"].get<std::string>();
	job.ff_dir = data["paths"]["ff_dir"].get<std::string>();
	
	return job;
}

// Parameters read_config(const std::string& filename){
// 	std::cout << "Start reading config file \n";
// 	std::ifstream file(filename);
// 	if (!file){
// 		throw std::runtime_error("Cannot open config file!");
// 	}

// 	std::unordered_map<std::string, std::string> values {};
// 	std::string line {};

// 	while (std::getline(file, line)){
// 		if (line.empty() || line[0] == '#'){
// 			continue;
// 		}
// 		std::istringstream iss(line);
// 		std::string key {};
// 		std::string eq {};
// 		std::string value {};

// 		iss >> key >> eq >> value;
// 		if (eq != "="){
// 			throw std::runtime_error("Invalid config format");
// 		}
// 		values[key] = value;
// 	}

// 	Parameters p {};
// 	p.N = std::stoi(values.at("N"));
// 	p.n_steps = std::stoi(values.at("n_steps"));
// 	p.dt = std::stod(values.at("dt")); //ps
// 	p.epsilon = std::stod(values.at("epsilon")); // kj/mol
// 	p.sigma = std::stod(values.at("sigma")); // nm
// 	p.mass_au = std::stod(values.at("mass")); // u
// 	p.pbc_L_nm = std::stod(values.at("pbc_L")); //nm, sidelength of simulation box centered around 0
// 	p.write_steps = std::stoi(values.at("write_steps"));
// 	p.output = values.at("output");
// 	p.output_opt = values.at("output_opt");
// 	p.ensemble_type = values.at("ensemble_type");
// 	p.init_pos = values.at("init_pos"); // way of initializing the positions
// 	p.init_dist = std::stod(values.at("init_dist"));
// 	p.T_K = std::stod(values.at("T_K"));
// 	p.optimization = std::stoi(values.at("optimization"));

// 	p.thermostat = parse_thermostat(values.at("thermostat"));

// 	std::cout << "---------------------------------------------------------------------------\n";
// 	std::cout << p.N << " atoms \n";
// 	std::cout << p.n_steps << " steps\n";
// 	std::cout << "timestep: " << p.dt << " ps \n";
// 	std::cout << "epsilon: " << p.epsilon << '\n';
// 	std::cout << "sigma: " << p.sigma << '\n';
// 	std::cout << "mass: " << p.mass_au << " g/mol \n";
// 	std::cout << "boxlength: " << p.pbc_L_nm << " nm \n";
// 	std::cout << "writing every: " << p.write_steps << " steps \n";
// 	std::cout << "ensemble type: " << p.ensemble_type << '\n';
// 	std::cout << "positions initialized by scheme: " << p.init_pos << '\n';
// 	std::cout << "initial distance between two particles: " << p.init_dist << " nm \n";
// 	std::cout << "T: " << p.T_K << " K \n";
// 	std::cout << "Optimization: " << p.optimization << '\n';
// 	std::cout << "Thermostat: " << values.at("thermostat") << '\n';
// 	std::cout << "Finished reading config \n"; 
// 	return p;
// }

Thermostat parse_thermostat(const std::string& key_thermostat){
	if (key_thermostat == "none"){
		return Thermostat::None;
	}
	if (key_thermostat == "berendsen"){
		return Thermostat::Berendsen;
	}
	if (key_thermostat == "vrescale"){
		return Thermostat::VRescale;
	}
	if (key_thermostat == "nosehoover"){
		return Thermostat::None;
	}
	throw std::runtime_error("Unknown thermostatL " + key_thermostat);
}

std::vector<std::string> get_forcefield_files(const std::string& ff_dir){
	std::cout << "Search directory " << ff_dir << " for forcefield definitions \n";
	std::vector<std::string> files {};
	for (const auto& entry: std::filesystem::directory_iterator(ff_dir)){
		if (entry.is_regular_file() && entry.path().extension() == ".json"){
			files.push_back(entry.path().string());
		}
	}
	std::cout << "Found forcefiled files:" << '\n';
	for (auto file: files){
		std::cout << file << '\n';
	}
	return files;
}


std::unordered_map<std::string, AtomType> parse_atom_types(const std::vector<std::string>& file_list){
	std::unordered_map<std::string, AtomType> types {};
	for (auto f: file_list){
		std::ifstream file(f);
		if (!file){
			throw std::runtime_error("Cannot open force field file");
		}
		json data;
		file >> data;
		for (auto [type_name, type_data]: data["types"].items()){
			AtomType atom {};
			atom.name = type_name;
			atom.charge = type_data["charge"];
			atom.element = type_data["element"];
			atom.mass = type_data["mass"];
			atom.sigma = type_data["sigma"];
			atom.epsilon = type_data["epsilon"];

			types.insert({atom.name, atom}); // already avoids double definition

		}
	}
	return types;
}

std::string make_bond_key(std::string string_a, std::string string_b){
	if (string_a < string_b){
		return string_a + ":" + string_b;
	}
	else {
		return string_b + ":" + string_a;
	}
}

std::vector<BondType> parse_bond_types(const std::vector<std::string>& file_list, const std::unordered_map<std::string, AtomType>& atom_types){
	std::vector<BondType> bond_types; 
	for (auto f: file_list){
		std::ifstream file(f);
		if (!file){
			throw std::runtime_error("Cannot open force field file");
		}
		json data;
		file >> data;
		for (auto bond_type: data["bonds"]){
			BondType bond {};
			bond.r0 = bond_type["r0"];
			bond.k = bond_type["k"];
			bond.name = make_bond_key(bond_type["atoms"].at(0), bond_type["atoms"].at(1));
			for (auto atom: bond_type["atoms"]){
				if (!atom_types.contains(atom)){
					throw std::runtime_error("Bonds make use of undefined atom type");
				}
			}
			if (!bondtype_exists(bond_types, bond.name)){
				bond_types.push_back(bond);
			}
		}
	}
	return bond_types;
}

bool bondtype_exists(std::vector<BondType> bond_types, std::string name_of_type){
	for (auto b: bond_types){
		if (b.name == name_of_type){
			return true;
		}
	}
	return false;
}

int get_type_idx(std::vector<BondType> bond_types, std::string name_of_type){
	for (size_t i=0; i<bond_types.size(); ++i){
		if (bond_types[i].name == name_of_type){
			return static_cast<int>(i);
		}
	}
	return -1;
}