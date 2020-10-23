//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....//
//                                                                     //
// Jet Physics code                                                    //
//                                                                     //
// Author: Omar Jamil                                                  //
// Contact: omar.jamil@gmail.com                                       //
// Astrophysical Institute                                             //
// Ohio University                                                     //
//                                                                     //
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....//



// The following file links the parameters and makes them available
// to all other classes. The variables declared here need to declared
// with extern keyword in other classes. The following are the
// definitions.

#include "parameters.hpp"

Parameters params("params.par");

int ext_x_cells = params.x_cells();
int ext_y_cells = params.y_cells();
int ext_z_cells = params.z_cells();
double ext_cell_size = params.cell_size();
int ext_dist_grid = params.dist_grid();
int ext_angle_dist_grid = params.angle_dist_grid();
double ext_e_gamma_min = params.e_gamma_Min();
double ext_e_gamma_max = params.e_gamma_Max();
double ext_e_gamma_brk = params.e_gamma_Brk();
double ext_pLaw_index = params.pLaw_index();
double ext_pLaw_index1 = params.pLaw_index1();
double ext_pLaw_norm = params.pLaw_norm();
double ext_kT_e = params.kT_e();
double ext_therm_norm = params.therm_norm();
double ext_ph_nu_min = params.ph_nu_min();
double ext_ph_nu_max = params.ph_nu_max();
std::string ext_B_file = params.b_file();
int ext_t_steps = params.t_steps();
bool ext_res_wrt_B = params.results_wrt_B();
double ext_source_d = params.source_dist();
double ext_source_z = params.source_redshift();
double ext_source_angle = params.s_inclination();
double ext_lorentz_factor = params.lorentz_factor();
int ext_debug_int = params.debug_int();
bool ext_compton_bool = params.compton_bool();
double ext_elec_ke_dens = params.elec_ke();
bool ext_use_norm_bool = params.use_norm();

