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
	initialize_random_positions(system, params);
	//initialize_positions_sc(system, params);
	// initialize_rand_velocity(system);
	// initialize_constant_velocity(system);
	// energy_opt(system, params);	
	std::ofstream fout_traj(params.output);
	std::ofstream fout_energy("energy_force.txt");
	for (int i = 0; i < params.n_steps; ++i){
		double Ekin = kinetic_energy(system, params);
		double Epot = potential_energy(system, params);
		double time = i * params.dt ;// time in ps
		write_energy_force(fout_energy, time, Ekin, Epot);
		velocity_verlet(system, params, params.dt);
		if (i % params.write_steps == 0){
			write_xyz(fout_traj, system);
		}
	}
	fout_traj.close();
	fout_energy.close();
	std::cout << "Finished simulation" << std::endl;
	return 0;
}



