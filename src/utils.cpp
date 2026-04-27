#include <fstream>
#include <cstdlib>
#include <iostream>
#include "system.h"
#include "utils.h"

void write_frame(std::ofstream& fout, const System& system){
	fout << system.particles.size() << "\n";
	fout << "Lattice=\"" << system.box_size	<< "0.0 0.0 0.0 " << system.box_size << " 0.0 0.0 0.0 " << system.box_size << '\"' << '\n';
	for (const auto& p: system.particles){
		fout << "Ar" << " " << p.x<< " " << p.y << " " << p.z << "\n";
	}
}

void write_frame_diagnostics(std::ofstream& fout, double time, double Ekin, double Epot, double T_K){
	fout << time << " " << Ekin << " " << Epot << " " << Ekin + Epot << " " << T_K << '\n';
}
