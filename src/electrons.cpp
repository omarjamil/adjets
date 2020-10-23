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
#include <vector>
#include <iostream>
#include <fstream>

#include "electrons.hpp"
#include "cell.hpp"
#include "libs/physcon.hpp"





#define sqr(a) ((a)*(a))
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
Electrons::Electrons(const int &size, const int &asize, const double &min,
		     const double &max):Distribution(size, asize, min, max), 
					gridS(size), agridS(asize), gMin(min), 
					gMax(max)
{
  momentumEnergy();
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
//copy constructor
// Electrons::Electrons(const Electrons *p) : Distribution(p), 
// 					   pPtr(p->pPtr)
// {}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
Electrons::~Electrons()
{
  delete pPtr;
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Electrons::momentumEnergy()
{
  pPtr = new std::vector<double>;
  
  for (std::vector<double>::iterator gamma = xPtr->begin();
       gamma != xPtr->end(); ++gamma)
    {
      pPtr->push_back(sqrt(sqr(*gamma) - 1.0));
    }
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Electrons::pNorm(const double &plindex, const double &eneD,
			const double &gamMin, const double &gamMax)
{
  double kappa;
  //eneDens from J/m^3 to gamma
  eneDens = eneD/(physcon.mc2());
  
  if(eneDens == 0.0)
    {
      kappa = 0.0;
    }
  else if(gamMax == 1.)
    {
      kappa = 0.0;
    }
  else
    {
    
      if(plindex == 2.0)
        {
          kappa = 
            eneDens/((log(gamMax) - log(gamMin)) + ((1./(gamMax)) 
                                                      - (1./(gamMin))));
               
        }
      else
        {
          
	  //use with SI units calc. of synchrotron
          kappa = /*pow(physcon.mc2(), (plindex-1.)) **/
            eneDens/(((1./(2.-plindex))*(pow(gamMax,(2.-plindex)) -
                                       pow(gamMin,(2.-plindex)))) -
                   ((1./(1.-plindex))*(pow(gamMax,(1.-plindex)) - 
                                       pow(gamMin,(1.-plindex)))));
          
        }
    }
  pNormalization = kappa;
  return pNormalization;
  
  
}


//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Electrons::powerLaw(const double ne, const double &p, 
			 const double &gammaMin, const double &volume)
{
  double n0 = ne*(p-1.0);
    
  for (std::vector<double>::iterator n = nPtr->begin(), 
	 gamma = xPtr->begin(); gamma != xPtr->end(); ++n, ++gamma) 
    {
      // quicker than starting at begin()+index(gamma_min)?
      if (*gamma >= gammaMin) 
	{
	  *n += n0 * pow(*gamma, -p);
	}
    }    
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Electrons::powerLawEnDen(const double ne, const double &p, const double &p1, 
			      const double &gammaMin, const double &gammaMax,
			      const double &gammaBrk, const double &volume)
{
  if (gammaBrk <= gammaMin || gammaBrk >= gammaMax)
    {
      std::cerr<<"\n Error! gamma_break should be >= gamma_min and <= gamma_max"
	       <<"\n";
      exit(1);
    }
   
  double n0 = pNorm(p, ne, gammaMin, gammaBrk);
  double scale = pow(gammaBrk,p1-p);
  
  
  for (std::vector<double>::iterator n = nPtr->begin(), 
	 gamma = xPtr->begin(); gamma != xPtr->end(); ++n, ++gamma) 
    {
      if (*gamma >= gammaMin && *gamma < gammaBrk) 
      	{
      	  *n += n0 * pow(*gamma, -p);
      	}
      else if (*gamma >= gammaBrk && *gamma <= gammaMax)
      	{
      	  *n += (n0 * scale) * pow(*gamma,-p1);
	}
    } 
  //Print the distribution
  // for (auto x = xPtr->begin(), y = nPtr->begin()
  // 	 ; x != xPtr->end(); ++x, ++y)
  //   {
  //     std::cout<<*x<<"\t"<<*y<<"\n";
  //   }
  
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Electrons::powerLawAlpha()
{
  int g;
  int alpha;
  std::vector<double>::const_iterator n, ga;
    
  for (g = 0, n = nPtr->begin(), ga = xPtr->begin(); g < nPtr->size();
       ++g, ++n, ++ga) 
    {
      for (alpha = 0; alpha < nPtr->size(); ++alpha) 
	{
	  nPtrPhi[alpha+(g*gridS)] = *n/(nPtr->size());
	}
    }

}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Electrons::powerLawPhiTheta()
{
  int g, phi, theta, face;
  std::vector<double>::const_iterator n, ga;
  double div, val;
  div = (this->agridSize()) * (this->agridSize());
  int eDistId = 6;
  
  for (g = 0, n = nPtr->begin(), ga = xPtr->begin(); g < nPtr->size();
       ++g, ++n, ++ga) 
    {
      for (phi = 0; phi < this->agridSize(); ++phi) 
	{
	  for (theta = 0; theta < this->agridSize(); ++theta)
	    {
	      val = (*n/div);
	      this->writeArr3Trans(g, phi, theta, val, eDistId);
	    }
	}
    }
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Electrons::thermal(const double &ne, const double &kTe)
{
  double thetae = kTe * 1000.0 * physcon.e / physcon.mc2();
  double gammae = thetae + 1.0;
  double pe = sqrt(sqr(gammae) - 1.0);
  
  if (thetae < 0.011) 
    { // i.e. p_e < 0.15 (use limiting function) with norm 
      double norm = ne * 0.79 * pow(thetae, -3.0/2.0); // number manipulated
      for (std::vector<double>::iterator n = nPtr->begin(), 
	     gamma = xPtr->begin(), p = pPtr->begin(); n != nPtr->end(); 
	   ++n, ++gamma, ++p) 
	{
	  *n += norm * *gamma * *p * exp((1.0-(*gamma))/thetae);
	}
    } 
  else 
    { // i.e. p_e >= 0.15 (use bessel function)
      double bes;
      double num;
      const double x = 1.0/thetae, xnu = 2.0;
      numerical.bessik(1.0/thetae, 2.0, bes);
      double norm = ne / (bes * thetae);
      for (std::vector<double>::iterator n = nPtr->begin(), 
	     gamma = xPtr->begin(), p = pPtr->begin(); 
	   n != nPtr->end(); ++n, ++gamma, ++p) 
	{
	  num = norm * *gamma * *p * exp(-(*gamma) / thetae);
	  if (num >= 1.0)
	    {
	      *n += num;
	    }
	  else
	    {
	      *n += 0.0;
	    }
	}
    }  
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Electrons::outputGammaRange(std::string &filename)
{
  std::vector<double>::const_iterator gamma;
  std::ofstream gammaFile(filename.c_str());
  
  for (gamma = xPtr->begin(); gamma != xPtr->end(); ++gamma)
    {
      gammaFile<<*gamma<<"\n";
    }
  gammaFile.close();
      
}
