#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "energy_opt.h"
#include "system.h"
#include "utils.h"

// std::vector<double> harmonic_force(double dx, double dy, double dz, const Parameters& params)
// {
//     double r2 { dx*dx + dy*dy + dz*dz };
// 	double r { std::sqrt(r2) };
// 	double inv_r { 1 / r };
//     if (r2 < 1e-12) return {0.0, 0.0, 0.0}; // TODO: use some maximal value, not 0

// 	double f {-2 * params.epsilon * (r - params.sigma)};
// 	double dx_norm { dx * inv_r };
// 	double dy_norm { dy * inv_r };
// 	double dz_norm { dz * inv_r };
//     return {f * dx_norm, f * dy_norm, f * dz_norm};
// }


// double pair_potential_harmonic(double dr2, Parameters params){
// 	double r { std::sqrt(dr2) };
// 	return params.epsilon * std::pow(r - params.sigma, 2);
// }
