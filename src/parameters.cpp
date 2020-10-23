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

#include <ctime>

#include "parameters.hpp"
#include "libs/physcon.hpp"

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
Parameters::Parameters(std::string fn)
{
  try 
    {
      loadParametersFile(fn);
    } 
  catch (const std::exception &e) 
    {
      throw;
    }
    
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
Parameters::~Parameters()
{}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Parameters::loadParametersFile(std::string filename)
  throw (const std::exception &)
{
  //just to output date...
  std::time_t rawtime;
  struct tm * timeinfo;
  
  time(&rawtime);
  timeinfo = localtime ( &rawtime );
  std::cout<<asctime(timeinfo)<<std::endl;
  //.....
  std::ifstream parameterFile(filename.c_str());
  if (!parameterFile.good())
    throw std::runtime_error("Parameter file failed to load");

  std::string label, dummy, dummy2;
  std::cout<<"Reading the parameters file...."<<"\n";
  
  while(!parameterFile.eof())
    {
      parameterFile >> label;
      
      if(label == "x_cells")
        {
          parameterFile >> dummy >> x_cells_;
	}
      else if(label == "y_cells")
        {
          parameterFile >> dummy >> y_cells_;
	}
      else if(label == "z_cells")
        {
          parameterFile >> dummy >> z_cells_;
	}
      else if(label == "cell_size")
        {
          parameterFile >> dummy >> cell_size_;
	}
      else if(label == "dist_grid_number")
        {
          parameterFile >> dummy >> dist_grid_;
	}
      else if(label == "angle_dist_grid")
        {
          parameterFile >> dummy >> angle_dist_grid_;
	}
      else if(label == "e_gamma_min")
        {
          parameterFile >> dummy >> e_gamma_min_;
	}
      else if(label == "e_gamma_max")
        {
          parameterFile >> dummy >> e_gamma_max_;
	}
      else if(label == "e_gamma_brk")
        {
          parameterFile >> dummy >> e_gamma_brk_;
	}
      else if(label == "pLaw_index_0")
        {
          parameterFile >> dummy >> pLaw_index_;
	}
      else if(label == "pLaw_index_1")
        {
          parameterFile >> dummy >> pLaw_index1_;
	}
      else if(label == "use_norm")
        {
          parameterFile >> dummy >> dummy2;
	  if(dummy2 == "y")
	    {
	      use_norm_bool_ = 1;
            }
          else
            {
              use_norm_bool_ = 0;
            }
	}
      else if(label == "pLaw_norm")
        {
          parameterFile >> dummy >> pLaw_norm_;
	}
      else if(label == "electron_ke_density")
        {
          parameterFile >> dummy >> elec_ke_;
	}
      else if(label == "kT_e")
        {
          parameterFile >> dummy >> kT_e_;
	}
      else if(label == "therm_norm")
        {
          parameterFile >> dummy >> therm_norm_;
	}
      else if(label == "ph_nu_min")
        {
          parameterFile >> dummy >> ph_nu_min_;
	}
      else if(label == "ph_nu_max")
        {
          parameterFile >> dummy >> ph_nu_max_;
	}
      else if(label == "B_field_file")
        {
          parameterFile >> dummy >> b_file_;
	}
      else if(label == "t_steps")
	{
	  parameterFile >> dummy >> t_steps_;
	}
       else if(label == "results_wrt_B")
        {
          parameterFile >> dummy >> dummy2;
	  if(dummy2 == "y")
	    {
	      results_wrt_B_ = 1;
            }
          else
            {
              results_wrt_B_ = 0;
            }
	}
      else if(label == "source_distance")
	{
	  parameterFile >> dummy >> distance_;
	}
      else if(label == "source_redshift")
	{
	  parameterFile >> dummy >> redshift_;
	}
      else if(label == "source_inclination")
	{
	  parameterFile >> dummy >> s_inclination_;
	}
      else if(label == "lorentz_factor")
	{
	  parameterFile >> dummy >> lorentz_;
	}
      else if(label == "debug_int")
	{
	  parameterFile >> dummy >> debug_int_;
	}
      else if(label == "compton_on")
	{
	  parameterFile >> dummy >> dummy2;
	  
	  if(dummy2 == "y")
	    {
	      compton_bool_ = 1;
            }
          else
            {
              compton_bool_ = 0;
            }
	}
      
    }
  parameterFile.close();
  
}

  
