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


//the run file containing the main function.

#include <iostream>
#include <iterator>
#include <ctime>
#include "simulation.hpp"
/*! \mainpage Modelling Blazar jets: Code documentation 
*/

/*! Blazar jet code.
     Currently in development stage, the following code aims to model
blazar jets. The code shall take the output from acceleration mechanisms 
simulations and use them as input for the lepton, and eventually hadron, 
populations. The radiative processes to be included are synchrotron, 
(self-)Compton scattering, and hadronic processes. Angle dependent synchrotron
has been implemented */
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

int main()
{
  
 
  std::time_t timeStart, timeEnd;
  time(&timeStart);

  Simulation run;
  run.create();
  
  time(&timeEnd);
  double diff = difftime(timeEnd, timeStart) / 60.;
  std::cout<<"CPU time usage for the Simualtion [minutes]: "<<diff<<std::endl;
  return 0;
  
}
