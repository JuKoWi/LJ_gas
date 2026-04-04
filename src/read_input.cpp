#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include "read_input.h"

Parameters read_config(const std::string& filename){
	std::cout << "Start reading config file" << std::endl;
	std::ifstream file(filename);
	if (!file){
		throw std::runtime_error("Cannot open config file!");
	}

	std::unordered_map<std::string, std::string> values;
	std::string line;

	while (std::getline(file, line)){
		if (line.empty() || line[0] == '#'){
			continue;
		}
		std::istringstream iss(line);
		std::string key, eq, value;

		iss >> key >> eq >> value;
		if (eq != "="){
			throw std::runtime_error("Invalid config format");
		}
		values[key] = value;
	}

	Parameters p;
	p.N = std::stoi(values.at("N"));
	p.n_steps = std::stoi(values.at("n_steps"));
	p.dt = std::stod(values.at("dt")); //ps
	p.epsilon = std::stod(values.at("epsilon")); // kj/mol
	p.sigma = std::stod(values.at("sigma")); // nm
	p.mass_au = std::stod(values.at("mass")); // u
	p.pbc_L_nm = std::stod(values.at("pbc_L")); //nm
	p.write_steps = std::stoi(values.at("write_steps"));
	p.output = values.at("output");
	p.output_opt = values.at("output_opt");

	std::cout << p.N << std::endl;
	std::cout << p.dt << std::endl;
	std::cout << p.n_steps << std::endl;

	std::cout << "Finished reading config" << std::endl;
	return p;
}


