#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include "read_input.h"
#include "system.h"

void initialize_random_positions(std::vector<Particle>& system, Parameters params){
	for (auto& p : system){
		p.x = drand48() *  params.pbc_L_nm - 0.5 * params.pbc_L_nm;
		p.y = drand48() *  params.pbc_L_nm - 0.5 * params.pbc_L_nm;
		p.z = drand48() *  params.pbc_L_nm - 0.5 * params.pbc_L_nm;
	}
}

void initialize_positions_sc(std::vector<Particle>& system, Parameters params){
	std::cout << "Starting position initializion" << std::endl;
	int N_direction = std::ceil(std::cbrt(params.N));
	double dr = params.pbc_L_nm / N_direction;
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

void initialize_constant_velocity(std::vector<Particle>& system){
	std::cout << "Initialize constant velocity" << std::endl;
	for (auto& p : system){
		p.vx = 1;
	}
	return;
}

void initialize_rand_velocity(std::vector<Particle>& system){
	std::cout << "Initialize random velocity" << std::endl;
	float vrange = 1;
	for (auto& p : system){
		p.vx =drand48() *  vrange - 0.5 * vrange;
	}
	return;
}

void initialize_positions_file(std::vector<Particle>& system){
	std::cout << "Initialize positions from file" << std::endl;
}


