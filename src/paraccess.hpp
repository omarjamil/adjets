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


#ifndef PARACCESS_HH
#define PARACCESS_HH

//enumerate the direction with integers
//first one is zero...
//the directions are wrt to the source i.e. away from 
//BH would be FORWARD
enum Direction {FORWARD, BACKWARD, UP, DOWN, RIGHT, LEFT};

//extern 'promises' the compiler that this will become available later
//including this file makes the following variable avaialable.
//this works in conjunction with parameters.* and paracess.cc

extern const int ext_x_cells;
extern const int ext_y_cells;
extern const int ext_z_cells;
extern const double ext_cell_size;
extern const int ext_dist_grid;
extern const int ext_angle_dist_grid;
extern const double ext_e_gamma_min;
extern const double ext_e_gamma_max;
extern const double ext_e_gamma_brk;
extern const double ext_pLaw_index;
extern const double ext_pLaw_index1;
extern const double ext_pLaw_norm;
extern const double ext_kT_e;
extern const double ext_therm_norm;
extern const double ext_ph_nu_min;
extern const double ext_ph_nu_max;
extern const std::string ext_B_file;
extern const int ext_t_steps;
extern const bool ext_res_wrt_B;
extern const double ext_source_d;
extern const double ext_source_z;
extern const double ext_source_angle;
extern const double ext_lorentz_factor;
extern const int ext_debug_int;
extern const bool ext_compton_bool;
extern const double ext_elec_ke_dens;
extern const bool ext_use_norm_bool;

#endif // PARACCESS_HH
