#ifndef PARAMETERS_H
#define PARAMETERS_H

enum class Thermostat{
	None,
	VRescale,
	Berendsen,
	NoseHoover,
};

struct Parameters{
	int N;
	int n_steps;
	double dt;
	double epsilon;
	double sigma;
	double mass_au;
	double pbc_L_nm;
	double T_K;
	int optimization;
	int write_steps;
	std::string output;
	std::string output_opt;
	std::string init_pos;
	double init_dist; 
	std::string ensemble_type;
	Thermostat thermostat;
};


Parameters read_config(const std::string& filename);
Thermostat parse_thermostat(const std::string& key_thermostat);

#endif
