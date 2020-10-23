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

#ifndef PARAMETERS_HH
#define PARAMETERS_HH

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
//! For reading the parameters file
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
class Parameters
{
  
public:
  Parameters(std::string); 	//the default constructor
  //Parameters(const Parameters &c);
  ~Parameters();
  
  void loadParametersFile(std::string) throw (const std::exception &);

  int x_cells()
  {
    return x_cells_;
  }
  int y_cells()
  {
    return y_cells_;
  }
  int z_cells()
  {
    return z_cells_;
  }
  double cell_size()
  {
    return cell_size_;
  }
  int dist_grid()
  {
    return dist_grid_;
  }
  
  int angle_dist_grid()
  {
    return angle_dist_grid_;
  }
  double e_gamma_Min()
  {
    return e_gamma_min_;
  }
  double e_gamma_Max()
  {
    return e_gamma_max_;
  }
  double e_gamma_Brk()
  {
    return e_gamma_brk_;
  }
  double pLaw_index()
  {
    return pLaw_index_;
  }
  double pLaw_index1()
  {
    return pLaw_index1_;
  }
  double pLaw_norm()
  {
    return pLaw_norm_;
  }
  double elec_ke()
  {
    return elec_ke_;
  }
  double kT_e()
  {
    return kT_e_;
  }
  double therm_norm()
  {
    return therm_norm_;
  }
  double ph_nu_min()
  {
    return ph_nu_min_;
  }
  double ph_nu_max()
  {
    return ph_nu_max_;
  }
  std::string b_file()
  {
    return b_file_;
  }
  
  int t_steps()
  {
    return t_steps_;
  }
  
  bool results_wrt_B()
  {
    return results_wrt_B_;
  }
  
  double source_dist()
  {
    return distance_;
  }
  
  double source_redshift()
  {
    return redshift_;
  }
  
  double s_inclination()
  {
    return s_inclination_;
  }

  double lorentz_factor()
  {
    return lorentz_;
  }
  
  //used for debugging
  int debug_int()
  {
    return debug_int_;
  }
  
  bool compton_bool()
  {
    return compton_bool_;
  }
  
  bool use_norm()
  {
    return use_norm_bool_;
  }
  
  
  
private:
  int x_cells_;
  int y_cells_;
  int z_cells_;
  double cell_size_;
  int dist_grid_;
  int angle_dist_grid_;
  double e_gamma_min_;
  double e_gamma_max_;
  double e_gamma_brk_;
  double pLaw_index_;
  double pLaw_index1_;
  bool use_norm_bool_;
  double pLaw_norm_;
  double elec_ke_;
  double kT_e_;
  double therm_norm_;
  double ph_nu_min_;
  double ph_nu_max_;
  std::string b_file_;
  int t_steps_;
  bool results_wrt_B_;
  double distance_;
  double redshift_;
  double s_inclination_;
  double lorentz_;
  int debug_int_;
  bool compton_bool_;
  
  

  

};
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

#endif // PARAMETERS_HH
