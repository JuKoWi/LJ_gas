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
	std::vector<BondType> bond_types = parse_bond_types(forcefield_files, atom_types);
	std::vector<AngleType> angle_types = parse_angle_types(forcefield_files, atom_types);

	for (auto [key, type]: atom_types){
		type.print();
	}

	for (auto type: angle_types){
		type.print();
	}

	System system {};
	system.initialize(job_params, atom_types, bond_types, angle_types);

	auto start_md { std::chrono::high_resolution_clock::now() };
	system.run_simulation(job_params);
	auto end_md {std::chrono::high_resolution_clock::now() };
	// fout_traj.close();
	// fout_energy.close();
	
	auto end_total { std::chrono::high_resolution_clock::now() };
	std::cout << "Finished simulation after " << std::chrono::duration<double>(end_total - start_total).count() << " s.\n";
	std::cout << "MD took " << std::chrono::duration<double>(end_md- start_md).count() << " s.\n";
	std::cout << "------------------------------------------------------------------------------------------------------------\n";  
	return 0;
}



