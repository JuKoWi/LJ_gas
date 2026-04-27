#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include "read_input.h"

Parameters read_config(const std::string& filename){
	std::cout << "Start reading config file \n";
	std::ifstream file(filename);
	if (!file){
		throw std::runtime_error("Cannot open config file!");
	}

	std::unordered_map<std::string, std::string> values {};
	std::string line {};

	while (std::getline(file, line)){
		if (line.empty() || line[0] == '#'){
			continue;
		}
		std::istringstream iss(line);
		std::string key {};
		std::string eq {};
		std::string value {};

		iss >> key >> eq >> value;
		if (eq != "="){
			throw std::runtime_error("Invalid config format");
		}
		values[key] = value;
	}

	Parameters p {};
	p.N = std::stoi(values.at("N"));
	p.n_steps = std::stoi(values.at("n_steps"));
	p.dt = std::stod(values.at("dt")); //ps
	p.epsilon = std::stod(values.at("epsilon")); // kj/mol
	p.sigma = std::stod(values.at("sigma")); // nm
	p.mass_au = std::stod(values.at("mass")); // u
	p.pbc_L_nm = std::stod(values.at("pbc_L")); //nm, sidelength of simulation box centered around 0
	p.write_steps = std::stoi(values.at("write_steps"));
	p.output = values.at("output");
	p.output_opt = values.at("output_opt");
	p.ensemble_type = values.at("ensemble_type");
	p.init_pos = values.at("init_pos"); // way of initializing the positions
	p.init_dist = std::stod(values.at("init_dist"));
	p.T_K = std::stod(values.at("T_K"));
	p.optimization = std::stoi(values.at("optimization"));

	p.thermostat = parse_thermostat(values.at("thermostat"));

	std::cout << "---------------------------------------------------------------------------\n";
	std::cout << p.N << " atoms \n";
	std::cout << p.n_steps << " steps\n";
	std::cout << "timestep: " << p.dt << " ps \n";
	std::cout << "epsilon: " << p.epsilon << '\n';
	std::cout << "sigma: " << p.sigma << '\n';
	std::cout << "mass: " << p.mass_au << " g/mol \n";
	std::cout << "boxlength: " << p.pbc_L_nm << " nm \n";
	std::cout << "writing every: " << p.write_steps << " steps \n";
	std::cout << "ensemble type: " << p.ensemble_type << '\n';
	std::cout << "positions initialized by scheme: " << p.init_pos << '\n';
	std::cout << "initial distance between two particles: " << p.init_dist << " nm \n";
	std::cout << "T: " << p.T_K << " K \n";
	std::cout << "Optimization: " << p.optimization << '\n';
	std::cout << "Thermostat: " << values.at("thermostat") << '\n';
	std::cout << "Finished reading config \n"; 
	return p;
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
