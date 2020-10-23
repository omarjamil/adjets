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

#ifndef SYNCHROTRON_HH
#define SYNCHROTRON_HH

#include "electrons.hpp"
#include "photons.hpp"
#include "libs/physcon.hpp"

//! Synch radiation class
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
class Synchrotron
{
  
public:
  
  Synchrotron(); 	//the default constructor
  Synchrotron(const Synchrotron &c);
  ~Synchrotron();
 
  //Pitch angle dependent Synchrotron
  void synchPhiTheta(Electrons *, Photons *, const double &, 
		     std::vector<double> &, 
		     const double &, const Arr4D<int, double> *);
   //Routine for calculate Modified bessel function
  double synRad(double &);

  double synRadAsymptote(double &);
  
protected:
  Numerical numerical;
 
  
private:
  
  
};
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....





#endif // SYNCHROTRON_HH
