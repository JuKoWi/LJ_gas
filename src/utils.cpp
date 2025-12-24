#include <fstream>
#include <cstdlib>
#include <iostream>
#include "system.h"

void write_xyz(std::ofstream& fout, const std::vector<Particle>& system){
	fout << system.size() << "\n\n";
	for (const auto& p: system){
		fout << "Ar" << " " << p.x << " " << p.y << " " << p.z << "\n";
	}
}
