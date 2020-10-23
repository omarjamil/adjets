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


#ifndef ELECTRONS_HH
#define ELECTRONS_HH

#include <vector>
#include "distribution.hpp"
#include "libs/numerical.hpp"
#include "libs/funcobj.hpp"

//! Electron distribution class
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
class Electrons : public Distribution
{
   
public:
  Electrons(const int &, const int &, const double &, const double &); 	//the default constructor
  //Electrons(const Electrons &c);
  //Electrons(const Electrons *);
  ~Electrons();
  //create a thermal electron distribution
  void thermal(const double &, const double &);
  //create a power-law electron distribution
  void powerLaw(const double, const double &, const double &, const double &);
  //power law using electron energyy density instead of normalization
  void powerLawEnDen(const double, const double &, const double &, 
		     const double &, const double &, const double &,
		     const double &);
  //create a 2D electron distribution which includes pitch angle
  void powerLawAlpha();
  //create a 3D electron distribution wrt pitch angle, and azimuthal
  void powerLawPhiTheta();
  
  void outputGammaRange(std::string &);
  
  inline double getPNorm()
  {
    return pNormalization;
  }

protected:
  Numerical numerical; //Electrons has-a Numerical object.
  void momentumEnergy();
  //calculate power law normalization from energy density
  double pNorm(const double &, const double &, const double &,
	       const double &);
  
  //returns electron velocity beta
  inline double beta(const double &gamma)
  {
    return sqrt(1. - (1. / pow(gamma, 2)));
  }
  

  
private:
  std::vector<double> *pPtr;  
  int gridS;
  int agridS;
  double gMin, gMax, eneDens;
  double pNormalization;
  GetIndex<double, double> getind;
};
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

#endif // ELECTRONS_HH



