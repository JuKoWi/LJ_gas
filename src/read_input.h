#ifndef PARAMETERS 
#define PARAMETERS

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

Parameters read_config(const std::string& filename);

#endif
