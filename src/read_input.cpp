#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>

struct Parameters{
	int N;
	int n_steps;
	double dt;
	double epsilon;
	double sigma;
	double mass_au;
	double pbcx_angstrom;
	double pbcy_angstrom;
	double pbcz_angstrom;
	int write_steps;
	std::string output;
};

Parameters read_config(const std::string& filename){
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
	p.dt = std::stod(values.at("dt"));
	p.epsilon = std::stod(values.at("epsilon"));
	p.sigma = std::stod(values.at("sigma"));
	p.mass_au = std::stod(values.at("mass"));
	p.pbcx_angstrom = std::stod(values.at("pbcx"));
	p.pbcy_angstrom = std::stod(values.at("pbcy"));
	p.pbcz_angstrom = std::stod(values.at("pbcz"));
	p.write_steps = std::stoi(values.at("write_steps"));
	p.output = values.at("output");

	return p;
}


