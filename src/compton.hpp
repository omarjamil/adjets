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

#ifndef COMPTON_HH
#define COMPTON_HH

#include "distribution.hpp"
#include "electrons.hpp"
#include "photons.hpp"

#define sqr(a) ((a)*(a))
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
class Compton
{
  
public:
  Compton(); 	//the default constructor
  Compton(const Compton &c);
  ~Compton();
  
  void comptonEmission(Electrons *, Photons *, double &);
 
  
protected:
  
private:
  double comptonKernel(double &, double &, double &, double &);
   
  Numerical *nums;

  //double comptoncsHeadOn(double &, double &, double &);
};
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....


#endif // COMPTON_HH
