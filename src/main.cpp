#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <string>
#include "read_input.h"
#include "system.h"
#include "utils.h"
#include "energy_opt.h"

int main(){
	std::cout << "Starting main" << std::endl;
	Parameters params = read_config("config.txt");
	std::vector<Particle> system(params.N);
	std::cout << "Initialized System" << std::endl;
	initialize_positions_sc(system, params);
	energy_opt(system, params);	
	std::ofstream fout(params.output);
	for (int i = 0; i < params.n_steps; ++i){
		velocity_verlet(system, params, params.dt);
		if (i % params.write_steps == 0){
			write_xyz(fout, system);
		}
	}
	fout.close();
	std::cout << "Finished simulation" << std::endl;
	return 0;
}



