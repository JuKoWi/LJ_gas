#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <chrono>
#include <string>
#include "read_input.h"
#include "system.h"
#include "utils.h"
#include "energy_opt.h"

int main(){
	auto start_total { std::chrono::high_resolution_clock::now() };

	std::cout << "Starting main" << std::endl;
	Parameters params = read_config("config.txt"); // kept constant except for box size
	System system {};
	system.initialize(params);
	std::cout << "Initialized System" << std::endl;
	if (params.optimization == 1){
		energy_opt(system, params);
	}
	std::ofstream fout_traj(params.output);
	std::ofstream fout_energy("energy_force.txt");
	double time {};
	double Ekin {};
	double Epot {};
	auto start_md { std::chrono::high_resolution_clock::now() };
	for (int i = 0; i < params.n_steps; ++i){
		velocity_verlet(system, params, params.dt);
		if (params.ensemble_type == "NVT"){
			system.apply_vrescale(params.T_K);
		}
		if (i % params.write_steps == 0){
			write_frame(fout_traj, system);
			time = i * params.dt ;// time in ps
			Ekin = kinetic_energy(system);
			Epot = potential_energy(system, params);
			write_frame_diagnostics(fout_energy, time, Ekin, Epot, system.get_temp());
		}
	}
	auto end_md {std::chrono::high_resolution_clock::now() };
	fout_traj.close();
	fout_energy.close();
	
	auto end_total { std::chrono::high_resolution_clock::now() };
	std::cout << "Finished simulation after " << std::chrono::duration<double>(end_total - start_total).count() << " s.\n";
	std::cout << "MD took " << std::chrono::duration<double>(end_md- start_md).count() << " s.\n";
	std::cout << "------------------------------------------------------------------------------------------------------------\n";  
	return 0;
}



