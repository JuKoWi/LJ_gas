#ifndef PARAMETERS_H
#define PARAMETERS_H

struct Parameters{
	int N;
	int n_steps;
	double dt;
	double epsilon;
	double sigma;
	double mass_au;
	double pbc_L_angstrom;
	int write_steps;
	std::string output;
	std::string output_opt;
};

Parameters read_config(const std::string& filename);

#endif
