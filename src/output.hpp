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

#ifndef OUTPUT_HH
#define OUTPUT_HH

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iterator>

#include "photons.hpp"
#include "electrons.hpp"
#include "cell.hpp"

//! Writing output files
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
class Output
{
public:
  Output();
  ~Output();
  //writes the photon files
  void writePhotonFile(Photons *, std::string &);
  //writes the electron files
  void writeElectronFile(Electrons *, std::string &);
  void writePhotElecOutput(std::vector<Cell *> *);
  void writeVolumeResults(Arr3D<int, double> *);
  void writeResArray(Arr3D<int, double> *, std::vector<Cell *> *, 
		     double &, double &, 
		     double &,  double &, int &);
  void outputArray(std::string &, std::vector<Cell *> *, 
		   Arr3D<int, double> *, double &, double &, 
		   double &, double &);
  void outputArrayAverage(std::string &, std::vector<Cell *> *, 
			  Arr3D<int, double> *, double &, 
			  double &, double &, double &);
  double bulkDopplerFactor(double &, double &);
  double dopplerFlux(double &, double &, double &);
  
  
};

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

#endif // OUTPUT_HH
