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

#ifndef PHOTONS_HH
#define PHOTONS_HH

#include <vector>
#include "distribution.hpp"
//! Photon class code
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
class Photons : public Distribution
{
  
public:
  Photons(const int &, const int &, const double &, const double &); 	//the default constructor
  //Photons(const Photons &c);
  //Photons(const Photons *p);
  ~Photons();

  void outputNuRange(std::string &, double &, double &);
  
protected:
  double nPhotons() const; //total number of photons
  
private:
  int gridS;
  int agridS;
  void initializeArr();
  
};
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

#endif // PHOTONS_HH
