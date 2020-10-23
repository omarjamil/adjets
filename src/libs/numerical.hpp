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

#ifndef NUMERICAL_HH
#define NUMERICAL_HH

#include <vector>
#include <cmath>
//! Various numerical techniques

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
class Numerical
{
  
public:
  Numerical(); 	//the default constructor
  //Numerical(const Numerical &c);
  ~Numerical();

  double chebev(const double, const double, std::vector<double> &, 
		const int, const double);
  void beschb(double, double *, double *, double *, double *);
  // void bessik(const double, const double, double &, double &, 
  // 	      double &, double &);
  //Modified Bessel function
  void bessik(const double, const double, double &);
  double gammalnFunction(const double &);
  double random(int &idum);
  //convert from sphreical to Cartesian coordinates
  void spheriToCarte(double &, double &, double &,
		     double &, double &, double &);
  //convert from Cartesian to spherical coordinates
  void carteToSpheri(double &, double &, double &,
		     double &, double &, double &);
  int whichQuadrant(double &, double &, double &);
  double pitchAngle(double &, double &, double &, double &);
  
  
protected:
  
private:
  
};

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....


#endif // NUMERICAL_HH
