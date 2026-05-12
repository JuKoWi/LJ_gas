#ifndef UTILS_H
#define UTILS_H

#include <fstream>
#include <cstdlib>
#include <iostream>
#include "system.h"
#include "utils.h"

//vmd expects angstrom instead of nm
void System::write_xyz_frame(std::ofstream& fout){
	fout << N << "\n";
	fout << "Lattice=\"" << box_size[0]	<< " 0.0 0.0 0.0 " << box_size[1] << " 0.0 0.0 0.0 " << box_size[2] << '\"' << '\n';
	for (int i=0; i<N; i++){
		fout << atom_type_dict.at(atom_types[i]).element << " " << x[i]*10 << " " << y[i]*10 << " " << z[i]*10 << "\n";
	}
}

//vmd expects angstrom instead of nm
void System::write_lammpsdump_frame(std::ofstream& fout, int step_count){
	fout << "ITEM: TIMESTEP \n";
	fout << step_count << '\n';
	fout << "ITEM: NUMBER OF ATOMS \n";
	fout << N << '\n';
	fout << "ITEM: BOX BOUNDS pp pp pp \n";
	fout << -box_size[0]*10/2 << " " << box_size[0]*10/2 << '\n';
	fout << -box_size[1]*10/2 << " " << box_size[1]*10/2 << '\n';
	fout << -box_size[2]*10/2 << " " << box_size[2]*10/2 << '\n';
	fout << "ITEM: ATOMS id type x y z \n";
	for (int i=0; i<N; i++){
		fout << i+1 << " " << atom_type_dict[atom_types[i]].element << " " << x[i]*10 << " " << y[i]*10 << " " << z[i]*10 << '\n';
	}
}

void System::write_frame_diagnostics(std::ofstream& fout, double time, double Ekin, double Epot, double T_K){
	fout << time << " " << Ekin << " " << Epot << " " << Ekin + Epot << " " << T_K << '\n';
}

#endif