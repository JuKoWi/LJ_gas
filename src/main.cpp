#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <chrono>
#include <string>
#include <unordered_map>
#include "read_input.h"
#include "system.h"
#include "utils.h"
#include "energy_opt.h"

int main(){
	auto start_total { std::chrono::high_resolution_clock::now() };

	std::cout << "Starting main" << std::endl;
	JobParameters job_params = parse_job("config.json"); 
	std::vector<std::string> forcefield_files = get_forcefield_files(job_params.ff_dir);
	std::unordered_map<std::string, AtomType> atom_types = parse_atom_types(forcefield_files);
	for (auto [key, type]: atom_types){
		type.print();
	}

	System system {};
	system.initialize_atoms(job_params.infile_geom_json, atom_types);
	for (int i = 0; i<system.N; i++){
		std::cout << system.atom_types[i] << '\n';
	}

	for (int i = 0; i<system.N; i++){
		std::cout << system.charges[i] << '\n';
	}
	std::cout << system.N << '\n';
	std::cout << "Initialized System" << std::endl;
	// if (params.optimization == 1){
	// 	energy_opt(system, params);
	// }
	auto start_md { std::chrono::high_resolution_clock::now() };
	system.run_simulation(job_params);
	// for (int i = 0; i < params.n_steps; ++i){
	// 	velocity_verlet(system, params, params.dt);
	// 	if (params.ensemble_type == "NVT"){
	// 		system.apply_vrescale(params.T_K);
	// 	}
	// 	if (i % params.write_steps == 0){
	// 		// write_xyz_frame(fout_traj, system);
	// 		write_lammpsdump_frame(fout_traj, system, i);
	// 		time = i * params.dt ;// time in ps
	// 		Ekin = kinetic_energy(system);
	// 		Epot = potential_energy(system, params);
	// 		write_frame_diagnostics(fout_energy, time, Ekin, Epot, system.get_temp());
	// 	}
	// }
	auto end_md {std::chrono::high_resolution_clock::now() };
	// fout_traj.close();
	// fout_energy.close();
	
	auto end_total { std::chrono::high_resolution_clock::now() };
	std::cout << "Finished simulation after " << std::chrono::duration<double>(end_total - start_total).count() << " s.\n";
	std::cout << "MD took " << std::chrono::duration<double>(end_md- start_md).count() << " s.\n";
	std::cout << "------------------------------------------------------------------------------------------------------------\n";  
	return 0;
}



