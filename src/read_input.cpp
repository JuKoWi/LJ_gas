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

std::string make_angle_key(std::string string_a, std::string string_b, std::string string_c){
	if (string_a < string_c){
		return string_a + ":" + string_b + ":" + string_c;
	}
	else {
		return string_c + ":" + string_b + ":" + string_a;
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
			if (!type_exists(bond_types, bond.name)){
				bond_types.push_back(bond);
			}
		}
	}
	return bond_types;
}

std::vector<AngleType> parse_angle_types(const std::vector<std::string>& file_list, const std::unordered_map<std::string, AtomType>& atom_types){
	std::vector<AngleType> angle_types;
	for (auto f: file_list){
		std::ifstream file(f);
		if (!file){
			throw std::runtime_error("Cannot open force field file");
		}
		json data;
		file >> data;
		for (auto angle_type: data["angles"]){
			AngleType angle {};
			angle.theta0_deg = angle_type["theta"];
			angle.k = angle_type["k"];
			angle.name = make_angle_key(angle_type["atoms"].at(0),angle_type["atoms"].at(1), angle_type["atoms"].at(2));
			for (auto atom: angle_type["atoms"]){
				if (!atom_types.contains(atom)){
					throw std::runtime_error("Angles make use of undefined atom type");
				}
			}
			if (!type_exists(angle_types, angle.name)){
				angle_types.push_back(angle);
			}

		}

	}
	return angle_types;
}
