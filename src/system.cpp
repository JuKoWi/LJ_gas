#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <random>
#include <cassert>
#include "json.hpp"
#include "read_input.h"
#include "system.h"
#include "utils.h"
#include "energy_opt.h"

using json = nlohmann::json;

void System::initialize(JobParameters job_params, std::unordered_map<std::string, AtomType> atom_types, std::vector<BondType> bond_types){
	initialize_atoms(job_params.infile_geom_json, atom_types);
	initialize_bonds(job_params.infile_geom_json, bond_types);
	energy_opt(job_params.outfile_opt_xyz);
    temperature_goal = job_params.T_K;
	initialize_temperature(job_params.T_K);
    update_forces();
}

void System::initialize_atoms(const std::string& geom_file, std::unordered_map<std::string, AtomType> type_dict){
    std::cout << "Initialize atom positions and connectivity \n";
    atom_type_dict = type_dict;
	std::ifstream file(geom_file);
	if (!file){
		throw std::runtime_error("Cannot open geometry file!");
	}
	json data;
	file >> data;

    for (auto L: data["box"]){
        box_size.push_back(L);
    }

    for (auto atom: data["atoms"]){
        atom_types.push_back(atom["type"]);
        x.push_back(atom["pos"][0]);
        y.push_back(atom["pos"][1]);
        z.push_back(atom["pos"][2]);
        x_at_verlet.push_back(atom["pos"][0]);
        y_at_verlet.push_back(atom["pos"][1]);
        z_at_verlet.push_back(atom["pos"][2]);
        vx.push_back(0.0);
        vy.push_back(0.0);
        vz.push_back(0.0);
        fx.push_back(0.0);
        fy.push_back(0.0);
        fz.push_back(0.0);
        sigma_nm.push_back(atom_type_dict.at(atom["type"]).sigma);
        charges.push_back(atom_type_dict.at(atom["type"]).charge);
        epsilon_kjmol.push_back(atom_type_dict.at(atom["type"]).epsilon);
        masses_gmol.push_back(atom_type_dict.at(atom["type"]).mass);
    }
    N = atom_types.size();
    verlet_list.resize(N);
    deg_of_freedom = 3 * N - 3;
    set_lj_pairs();

    // set lj and verlet cutoffs
    // lj based on maximal pair-sigma
    // constant verlet skin
    double max_pair_cutoff {};
    for (int i=0; i<N; ++i){
        for (int j=0; j<N; ++j){
            double cutoff = 2.5 * pair_sigma_nm[i*N + j];
            max_pair_cutoff = std::max(cutoff, max_pair_cutoff);
        }
    }
    lj_cutoff_nm = max_pair_cutoff;
    verlet_cutoff_nm = lj_cutoff_nm + verlet_skin_nm;
    update_verlet();

    assert(atom_types.size() == x.size());
    assert(x.size() == y.size());
    assert(y.size() == z.size());
    assert(x.size() == fx.size());
    assert(x.size() == fy.size());
    assert(x.size() == fz.size());
    assert(pair_epsilon_kjmol.size() == x.size()*x.size());
    assert(pair_sigma_nm.size() == x.size()*x.size());

	std::cout << "Initialized System" << std::endl;
}


void System::initialize_bonds(const std::string& geom_file, std::vector<BondType> types_from_file){
    bond_types = types_from_file;

	std::ifstream file(geom_file);
	if (!file){
		throw std::runtime_error("Cannot open geometry file!");
	}
	json data;
	file >> data;
    for (auto bond: data["bonds"]){
        int atom_idx1 = bond[0];
        int atom_idx2 = bond[1];
        std::string bond_name {make_bond_key(atom_types[atom_idx1], atom_types[atom_idx2])};
        if (!bondtype_exists(types_from_file, bond_name)){
            throw std::runtime_error("No bond parameters defined for this combination of atom types");
        }
        bonds.i.push_back(atom_idx1);
        bonds.j.push_back(atom_idx2);
        bonds.type.push_back(get_type_idx(types_from_file, bond_name));
    }
    
}

void System::set_lj_pairs(){
    pair_epsilon_kjmol.resize(N * N);
    pair_sigma_nm.resize(N * N);
    for (int i=0; i<N; ++i){
        for (int j=0; j<N; ++j){
            pair_sigma_nm[i*N + j] = 0.5 * (sigma_nm[i] + sigma_nm[j]);
            pair_epsilon_kjmol[i*N + j] = std::sqrt(epsilon_kjmol[i] * epsilon_kjmol[j]);
        }
    }
}

void System::do_verlet_step(double dt){
    for (int i=0; i<N; i++){
        vx[i] += 0.5 * fx[i] / masses_gmol[i] * dt;
        vy[i] += 0.5 * fy[i] / masses_gmol[i] * dt;
        vz[i] += 0.5 * fz[i] / masses_gmol[i] * dt;

        x[i] += vx[i] * dt;
        y[i] += vy[i] * dt;
        z[i] += vz[i] * dt;
    }
    update_positions();
    update_forces();
    for (int i=0; i<N; i++){
        vx[i] += 0.5 * fx[i] / masses_gmol[i] * dt;
        vy[i] += 0.5 * fy[i] / masses_gmol[i] * dt;
        vz[i] += 0.5 * fz[i] / masses_gmol[i] * dt;
    }
}

void System::do_simulation_step(double dt){ // step_count as additional argument
    do_verlet_step(dt);
    apply_vrescale(temperature_goal);
}

// distance components while considering PBC
double minimum_image_distance(double dx, double L){
	if (dx > 0.5 * L){
		dx -= L;
	}
	if (dx < -0.5 * L){
		dx += L;
	}
	return dx;
}

void System::update_forces(){
    std::fill(fx.begin(), fx.end(), 0.0);
    std::fill(fy.begin(), fy.end(), 0.0);
    std::fill(fz.begin(), fz.end(), 0.0);

    update_lj_forces(true);
    update_bond_forces();
}

void System::update_lj_forces(bool verlet){
	if (verlet){
        if (verlet_update_required()){
            update_verlet();
        }
		for (int i = 0; i<N; ++i){
			for (auto& j : verlet_list[i]){
				if (j<=i){
					continue;
				}
				double dx { minimum_image_distance(x[j] - x[i], box_size[0]) };
				double dy { minimum_image_distance(y[j] - y[i], box_size[1]) };
				double dz { minimum_image_distance(z[j] - z[i], box_size[2]) };

				if (dx*dx + dy*dy + dz*dz < lj_cutoff_nm * lj_cutoff_nm){
                    Vec3 f_lj {lj_force(dx, dy, dz, pair_sigma_nm[i*N + j], pair_epsilon_kjmol[i*N +j])};
				    fx[i] -= f_lj.x;
				    fy[i] -= f_lj.y;
				    fz[i] -= f_lj.z;

				    fx[j] += f_lj.x;
				    fy[j] += f_lj.y;
				    fz[j] += f_lj.z;
				}
			}
		}
	}
    else{
        for (int i=0; i<N; ++i){
            for (int j=i+1; j<N; ++j){
                double dx {minimum_image_distance(x[j] - x[i], box_size[0])};
                double dy {minimum_image_distance(y[j] - y[i], box_size[1])};
                double dz {minimum_image_distance(z[j] - z[i], box_size[2])};
                if (dx*dx + dy*dy + dz*dz < std::pow(lj_cutoff_nm, 2)){
                    Vec3 f_lj {lj_force(dx, dy, dz, pair_sigma_nm[i*N + j], pair_epsilon_kjmol[i*N +j])};
				    fx[i] -= f_lj.x;
				    fy[i] -= f_lj.y;
				    fz[i] -= f_lj.z;

				    fx[j] += f_lj.x;
				    fy[j] += f_lj.y;
				    fz[j] += f_lj.z;
                }
            }
        }
    }
}

void System::update_bond_forces(){
    for (size_t idx=0; idx<bonds.i.size(); ++idx){
        double dx {minimum_image_distance(x[bonds.j[idx]] - x[bonds.i[idx]], box_size[0])};
        double dy {minimum_image_distance(y[bonds.j[idx]] - y[bonds.i[idx]], box_size[1])};
        double dz {minimum_image_distance(z[bonds.j[idx]] - z[bonds.i[idx]], box_size[2])};
        Vec3 f_bond {bond_force(dx, dy, dz,
            bond_types[bonds.type[idx]].r0,
            bond_types[bonds.type[idx]].k)
        };
        fx[bonds.i[idx]] += f_bond.x;
        fy[bonds.i[idx]] += f_bond.y;
        fz[bonds.i[idx]] += f_bond.z;

        fx[bonds.j[idx]] -= f_bond.x;
        fy[bonds.j[idx]] -= f_bond.y;
        fz[bonds.j[idx]] -= f_bond.z;
    }
}

inline Vec3 bond_force(double dx, double dy, double dz, double r0, double k){
    Vec3 force {};
    double dr {std::sqrt(dx*dx + dy*dy + dz*dz)};
    double displacement {dr - r0};
    double abs_val_force {k * displacement};
    force.x = abs_val_force * dx/dr;
    force.y = abs_val_force * dy/dr;
    force.z = abs_val_force * dz/dr;
    return force;
}


// Return components of LJ-force between two particles
inline Vec3 lj_force(double dx, double dy, double dz, double sigma, double epsilon) {
    double r2 {dx*dx + dy*dy + dz*dz};
    Vec3 f_vec {};
    if (r2 < 1e-6){
        r2 = 1e-6;
    } 
	double r {std::sqrt(r2)};

	double inv_r {1 / r};
	double f {24 * epsilon * (2 * std::pow(sigma * inv_r, 12) * inv_r - std::pow(sigma * inv_r, 6) * inv_r)}; //abs value force
	double dx_norm {dx * inv_r};
	double dy_norm {dy * inv_r};
	double dz_norm {dz * inv_r};
    f_vec.x = f * dx_norm;
    f_vec.y = f * dy_norm; 
    f_vec.z = f * dz_norm;
    return f_vec;
}

// map positions outside simulation box back to simulation box
void System::update_positions(){ // 0-centered simulation box
    double half_box_x {box_size[0]/2};
    double half_box_y {box_size[1]/2};
    double half_box_z {box_size[2]/2};
	for (int i=0; i<N; ++i){
		x[i] -= box_size[0] * std::floor((x[i] + half_box_x)/ box_size[0]);
		y[i] -= box_size[1] * std::floor((y[i] + half_box_y)/ box_size[1]);
		z[i] -= box_size[2] * std::floor((z[i] + half_box_z)/ box_size[2]);
	}
}

void System::run_simulation(JobParameters params){
	std::ofstream fout_traj(params.outfile_traj_lammps);
    std::ofstream fout_energy(params.outfile_energy);
	for (int i = 0; i < params.n_steps; ++i){
		do_simulation_step(params.dt);
		if (i % params.write_steps == 0){
			write_lammpsdump_frame(fout_traj, i);
            write_frame_diagnostics(fout_energy, 
                i * params.dt, 
                get_kinetic_energy(), 
                get_total_potential(), 
                get_temp());
			// write_xyz_frame(fout_traj, system);
		}
	}
}

double System::get_kinetic_energy(){
    double Ekin {};
    for (int i=0; i<N; ++i){
		Ekin += masses_gmol[i]/2 * vx[i] * vx[i]; // component wise contribution in kj/mol
		Ekin += masses_gmol[i]/2 * vy[i] * vy[i]; // component wise contribution in kj/mol
		Ekin += masses_gmol[i]/2 * vz[i] * vz[i]; // component wise contribution in kj/mol
    }
    return Ekin;
}

// // LJ potential contribution for one pair in kj/mol
double System::lj_pair_potential(double dr_squared, int id1, int id2){
	double sqr_sigma { pair_sigma_nm[id1*N + id2] * pair_sigma_nm[id1*N + id2]};
	return pair_epsilon_kjmol[id1*N + id2]+ 4 * pair_epsilon_kjmol[id1*N + id2] * (std::pow(sqr_sigma / dr_squared, 6) - std::pow(sqr_sigma / dr_squared, 3));
}

double System::bond_pair_potential(double dr, double r0, double k){
    return k/2 * (dr - r0) * (dr - r0);
}

// total potential energy in kj/mol
double System::get_total_lj_potential(){
	double dr_squared {};
	double dx {};
	double dy {};
	double dz {};
	double U { 0 };
	for (int i=0; i<N; ++i){
		for (int j=i+1; j<N; ++j){
			dx = minimum_image_distance(x[j] - x[i], box_size[0]);
			dy = minimum_image_distance(y[j] - y[i], box_size[1]);
			dz = minimum_image_distance(z[j] - z[i], box_size[2]);
			dr_squared = dx * dx + dy * dy + dz * dz;
			U += lj_pair_potential(dr_squared, i, j);
		}
	}
	return U;
}

double System::get_total_bond_potential(){
    double U {};
    double dx {};
    double dy {};
    double dz {};
    double dr {};
    for (size_t idx=0; idx<bonds.i.size(); ++idx){
        dx = minimum_image_distance(x[bonds.j[idx]] - x[bonds.i[idx]], box_size[0]);
        dy = minimum_image_distance(y[bonds.j[idx]] - y[bonds.i[idx]], box_size[1]);
        dz = minimum_image_distance(z[bonds.j[idx]] - z[bonds.i[idx]], box_size[2]);
        dr = std::sqrt(dx*dx + dy*dy + dz*dz);
        U += bond_pair_potential(dr, bond_types[bonds.type[idx]].r0, bond_types[bonds.type[idx]].k);
    }
    return U;
}

double System::get_total_energy(){
    return get_kinetic_energy() + get_total_potential();
}

double System::get_total_potential(){
    return get_total_lj_potential() + get_total_bond_potential();
}

double System::get_temp(){
	return 2.0 * get_kinetic_energy() / (deg_of_freedom * PhysicalConstants::K_B);
}

double System::get_square_distance(int i, int j){
    double dx {minimum_image_distance(x[j] - x[i], box_size[0])};
    double dy {minimum_image_distance(y[j] - y[i], box_size[1])};
    double dz {minimum_image_distance(z[j] - z[i], box_size[2])};
    return dx * dx + dy * dy + dz * dz;
}

void System::update_verlet(){
	for (int i=0; i<N; i++){
		verlet_list[i].clear();
		for (int j = 0; j < N; j++){
			if (i == j){
				continue;
			}
			if (get_square_distance(i,j) < std::pow(verlet_cutoff_nm, 2)){
				verlet_list[i].push_back(j);
			}
		}
	}
	for (int i=0; i<N; i++){
                  x_at_verlet[i] = x[i];
                  y_at_verlet[i] = y[i];
                  z_at_verlet[i] = z[i];
	}
}

bool System::verlet_update_required(){
    double max_disp2 {0};
	for (int i=0; i<N; i++){
        double dx {minimum_image_distance(x_at_verlet[i] - x[i], box_size[0])};
        double dy {minimum_image_distance(y_at_verlet[i] - y[i], box_size[1])};
        double dz {minimum_image_distance(z_at_verlet[i] - z[i], box_size[2])};
        max_disp2 = std::max(max_disp2, dx*dx + dy*dy + dz*dz);
    }
    return (max_disp2 > 0.25 * verlet_skin_nm * verlet_skin_nm);
}

// steepest descend energy optimization
void System::energy_opt(std::string opt_out){
    std::cout << "Start energy optimization " << '\n';
    std::ofstream fout(opt_out);

    double max_epsilon_kjmol {0};
    for (int i=0; i<N; ++i){
        max_epsilon_kjmol = std::max(max_epsilon_kjmol, epsilon_kjmol[i]);
    }
	double threshold { 1e-7 * max_epsilon_kjmol * N };

	double total_force_sqr { 0 };
	double eta { 1 };
	double maxstep { 0.01 }; //nm
	int max_rep { 1000 };
	write_xyz_frame(fout);
	double Epot_before {};
	double Epot_after {};

	for (int i=0; i<max_rep; ++i){
		total_force_sqr = 0;
		update_forces();
        for (int i=0; i<N; ++i){
            total_force_sqr += fx[i] * fx[i];
            total_force_sqr += fy[i] * fy[i];
            total_force_sqr += fz[i] * fz[i];
        }
		if (i % 100 == 0){
			std::cout << "Iteration: " << i << " Total force: " << total_force_sqr << '\n';
		}
		if (total_force_sqr < threshold){
			std::cout << "Minimized energy below force-threshold \n";
			fout.close();
			return;
		}

		Epot_before = get_total_potential();
        for (int i=0; i<N; ++i){
            double dx {eta * fx[i]};
            x[i] += (std::abs(dx) < maxstep) ? dx : std::copysign(maxstep, dx);
            double dy {eta * fy[i]};
            y[i] += (std::abs(dy) < maxstep) ? dy : std::copysign(maxstep, dy);
            double dz {eta * fz[i]};
            z[i] += (std::abs(dz) < maxstep) ? dz : std::copysign(maxstep, dz);
        }
		update_positions();
		write_xyz_frame(fout);
		Epot_after = get_total_potential();
		if (Epot_after > Epot_before){
			eta *= 0.5;
		}
		else{
			eta *= 1.1;
		}
	}
	fout.close();
	std::cout << "Stopped energy minimization after 1000 iterations \n";
}

void System::initialize_temperature(double T_K){
	std::random_device rd;
	std::mt19937 gen(rd());

    for (int i=0; i<N; ++i){
        double std_deviation = std::sqrt(PhysicalConstants::K_B * T_K / masses_gmol[i]);
        std::normal_distribution<double> dist(0.0, std_deviation);
        vx[i] = dist(gen);
        vy[i] = dist(gen);
        vz[i] = dist(gen);
    }

    //remove drift
    Vec3 com_momentum {get_cm_momentum()};
    double total_mass {};
    for (int i=0; i<N; ++i){
        total_mass += masses_gmol[i];
    }
    double v_com_x = com_momentum.x/total_mass;
    double v_com_y = com_momentum.y/total_mass;
    double v_com_z = com_momentum.z/total_mass;
    for (int i=0; i<N; ++i){
        vx[i] -= v_com_x;
        vy[i] -= v_com_y;
        vz[i] -= v_com_z;
    }
    apply_vrescale(T_K);
}

void System::apply_vrescale(double T_K){
	if (T_K == 0){
        for (int i=0; i<N; ++i){
            vx[i] = 0.0;
            vy[i] = 0.0;
            vz[i] = 0.0;
        }
	}
	else{
		double scaling_factor { std::sqrt(T_K / get_temp()) };
        for (int i=0; i<N; ++i){
            vx[i] *= scaling_factor;
            vy[i] *= scaling_factor;
            vz[i] *= scaling_factor;
        }
	}
}

// get center-of-mass momentum in g * nm/(mol * ps)
Vec3 System::get_cm_momentum(){
    Vec3 com_momentum {};
	for (int i=0; i<N; ++i){
		com_momentum.x += masses_gmol[i] * vx[i];
		com_momentum.y += masses_gmol[i] * vy[i];
		com_momentum.z += masses_gmol[i] * vz[i];
	}
	return com_momentum;
}


// // double System::get_pressure(){
// // 	double V {std::pow(box_size, 3)};
// // }

// void System::initialize(const Parameters& params){
// 	box_size = params.pbc_L_nm;
// 	temperature = params.T_K;

// 	x.resize(params.N);
// 	y.resize(params.N);
// 	z.resize(params.N);

// 	x_at_verlet.resize(params.N);
// 	y_at_verlet.resize(params.N);
// 	z_at_verlet.resize(params.N);

// 	vx.resize(params.N);
// 	vy.resize(params.N);
// 	vz.resize(params.N);

// 	fx.resize(params.N);
// 	fy.resize(params.N);
// 	fz.resize(params.N);

// 	masses_gmol.resize(params.N);
// 	charges.resize(params.N);

// 	N = params.N;
// 	deg_of_freedom = 3 * N - 3 ;


// }
