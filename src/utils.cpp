#include <fstream>
#include <cstdlib>
#include <iostream>
#include "system.h"
#include "utils.h"

void write_xyz(std::ofstream& fout, const std::vector<Particle>& system){
	fout << system.size() << "\n\n";
	for (const auto& p: system){
		fout << "Ar" << " " << p.x << " " << p.y << " " << p.z << "\n";
	}
}

void write_energy_force(std::ofstream& fout, double time, double Ekin, double Epot){
	fout << time << " " << Ekin << " " << Epot << " " << Ekin + Epot << "\n";
}
