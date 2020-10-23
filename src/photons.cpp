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

#include <cmath>
#include <algorithm>
#include <fstream>

#include "photons.hpp"
#include "libs/physcon.hpp"


#define sqr(a) ((a)*(a))

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
Photons::Photons(const int &size, const int &asize, const double &min,
		 const double &max) : Distribution(size, asize, min, max), 
				      gridS(size), agridS(asize)
{
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
Photons::~Photons()
{
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Photons::initializeArr()
{
  int i, j, k;
  double val = 0.0;
  
  for (i = 0; i < gridS; ++i)
    {
       for (j = 0; j < agridS; ++j)
	 {
	   for (k = 0; k < agridS; ++k)
	     {
	       //std::cout<<i<<"\t"<<j<<"\t"<<k<<"\t"<<std::endl;
	       
	     }
	 }
    }
}


//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Photons::nPhotons() const
{
  double integral = 0.0, constN = 4.0 * physcon.pi / (physcon.h * physcon.c);
  
  for (std::vector<double>::const_iterator iNu = nPtr->begin(), 
	 nu = xPtr->begin(), dnu = dxPtr->begin(); iNu != nPtr->end(); 
       ++iNu, ++nu, ++dnu)
    {
      integral += *iNu / *nu * *dnu;
    }
  
  return integral*constN;
} 

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Photons::outputNuRange(std::string &filename, double &redshift, 
			    double &doppVal)
{
  std::vector<double>::const_iterator nu;
  std::ofstream nuFile(filename.c_str());
  double freqShift = doppVal/(1.+redshift);
  
  for (nu = xPtr->begin(); nu != xPtr->end(); ++nu)
    {
      nuFile<<*nu<<"\t"<<*nu*freqShift<<"\n";
    }
  nuFile.close();
      
}
