#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include "read_input.h"
#include "system.h"

void initialize_random_positions(std::vector<Particle>& system, Parameters params){
	for (auto& p : system){
		p.x = drand48() *  params.pbc_L_angstrom - 0.5 * params.pbc_L_angstrom;
		p.y = drand48() *  params.pbc_L_angstrom - 0.5 * params.pbc_L_angstrom;
		p.z = drand48() *  params.pbc_L_angstrom - 0.5 * params.pbc_L_angstrom;
	}
}

void initialize_positions_sc(std::vector<Particle>& system, Parameters params){
	std::cout << "Starting position initializion" << std::endl;
	int N_direction = std::ceil(std::cbrt(params.N));
	double dr = params.pbc_L_angstrom / N_direction;
	int count = 0;
	for (int i=0; i<N_direction; ++i){
		for (int j=0; j<N_direction; ++j){
			for (int k=0; k<N_direction; ++k){
				if (count == params.N){
					return;
				}
				system[count].x = i * dr;
				system[count].y = j * dr;
				system[count].z = k * dr;
				count += 1;
			}
		}
	}
}


