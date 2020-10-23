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
#include <vector>

#include "compton.hpp"
#include "libs/physcon.hpp"
#include "libs/numerical.hpp"

//#define sqr(a) ((a)*(a))


//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
Compton::Compton()
{
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
Compton::~Compton()
{
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Compton::comptonEmission(Electrons *ele, Photons *pho, double &length)
{
  double c = 2.99792458e8;
  double mc2 = 9.1093897e-31*sqr(c);
  double h = 6.6260755e-34;
  double pi = 3.141592653589793;
  double r_e = 2.8179402894e-15;
  double sigma_t = 6.652462e-29;
  
  //  std::cout<<"Compton Scattering...";
  int  itScatNu, itElePhi, itEleTheta, itPhoPhi, itPhoTheta, itEleGamma, 
    itPhoNu, face;
  int gridS = pho->gridSize();
  int agridS = pho->agridSize();
  int intrinsic = 6;
  
  double eleTheta, elePhi, eleGamma, eleBeta, phoTheta, phoPhi, phoNu, phoNuS, 
    phoEpsilon, phoEpsilonLimit, phoEpsilonS, ERF_factor, inPhotE_ERF, outPhotE_ERF, nEle, 
    nPhotNu, nPhotTheta, nPhotPhi, dGamma, dNu, iNu;
  double dopp, muPh, muEl,
    cosPsi = 0.0, cosChi = 0.0, cosChi_ERF=0.0, cosThetaERF=0.0, crossSection = 0.0, intrInu=0.0;
  double r = 1.0, jComp = 0.0, jComp2 = 0.0;
  double energyFactor = h/mc2;
  double dTheta = 2.*pi/(agridS-1);
  double dPhi = 2.*pi/(agridS-1);
  double dCosPhi=cos(dPhi), cosPhi;
  double yParam, lowerLimit, upperLimit, thompsonDepth;
  
   //loop over scattered photon frequency
  for (itScatNu = 0; itScatNu <= gridS-1; ++itScatNu)
    {
      //scattering photon frequency
      phoNuS = pho->xVal(itScatNu);
      //scattering photon energy in electron rest mass units
      phoEpsilonS = phoNuS*energyFactor;
      
      //loop over electron phi; this is also scattered photon phi
      for (itElePhi = 0; itElePhi <= agridS-1; ++itElePhi)
	{
	  elePhi = ele->phiElement(itElePhi);
	   muEl = cos(elePhi);
	  //loop over electron theta; this is also scattered photon theta
	  for (itEleTheta = 0; itEleTheta <= agridS-1; ++itEleTheta)
	    {
	      eleTheta = ele->thetaElement(itEleTheta);
	      //distribution to which the photons will be passed
	      face = pho->getDirection(itScatNu, itElePhi, itEleTheta);
	      
	      //loop over photon phi to integrate
	      for (itPhoPhi = 0; itPhoPhi <= agridS-1; ++itPhoPhi)
		{
		  phoPhi = pho->phiElement(itPhoPhi);
		  cosPhi=cos(phoPhi);
		  muPh = cos(phoPhi);		  
		  //loop over photon theta to integrate
		  for (itPhoTheta = 0; itPhoTheta <= agridS-1; ++itPhoTheta)
		    {
		      phoTheta = pho->thetaElement(itPhoTheta);
		      		      
		      cosPsi = muPh*muEl + sqrt(1-sqr(muPh))*sqrt(1-sqr(muEl))*cos(phoTheta- eleTheta);
		      //integrate over interacting photon frequency
		      for(itPhoNu = 0; itPhoNu <= gridS-1; ++itPhoNu)
			{
			  phoNu = pho->xVal(itPhoNu);
			  phoEpsilon = phoNu * energyFactor;//interacting photon energy in mc2
			  nPhotNu = pho->get3DdistElement(itPhoNu, itPhoPhi, itPhoTheta, intrinsic);
			  dNu = pho->dxVal(itPhoNu);
			  //integrate over electron energy
			  for (itEleGamma = 0; itEleGamma <= gridS-1; ++itEleGamma)
			    {
			      eleGamma = ele->xVal(itEleGamma);
			      eleBeta =  (pow(eleGamma, 2) - 1)/(pow(eleGamma, 2));
			      nEle = ele->get3DdistElement(itEleGamma, itElePhi, itEleTheta, intrinsic);
			      dGamma = ele->dxVal(itEleGamma);
			      dopp = 1.0 - (eleBeta*cosPsi);
			      ERF_factor = eleGamma * (1.0 - eleBeta*cosPsi);//e rest frame
					      
			      inPhotE_ERF = phoEpsilon * ERF_factor;
			      outPhotE_ERF = phoEpsilonS * eleGamma * (1.0 - eleBeta*muEl);
			      cosChi_ERF = 1.+ (1./inPhotE_ERF) - (1./outPhotE_ERF);
			      cosThetaERF = (muEl - eleBeta)/(1.-eleBeta*muEl);
			      yParam = 1. - (phoEpsilonS/eleGamma);
			      //limit to stop scattered photon having more energy than the \
			      //the scattering electron
			      cosChi = 1. + (1./inPhotE_ERF) - (1./outPhotE_ERF);
			      if(cosChi > 1.) cosChi = 1.;
			      phoEpsilonLimit = inPhotE_ERF/(1. + (inPhotE_ERF*(1.-cosChi)));
			      
			      //lowerLimit = 1./(1.+(2.*inPhotE_ERF));
			      lowerLimit = inPhotE_ERF/2./eleGamma;
			      upperLimit = (2.*eleGamma*inPhotE_ERF)/(1.+2.*inPhotE_ERF);
			     			      
			      if(phoEpsilonS > lowerLimit && phoEpsilonS < upperLimit)
				{
				  crossSection = (pi*sqr(r_e)/(eleGamma*inPhotE_ERF))
				    *comptonKernel(phoEpsilonS,inPhotE_ERF, yParam, eleGamma);
				}
			      else
				{
				  crossSection = 0.0;
				}
			      			      
			      jComp += crossSection * nEle * dGamma;
			      crossSection = 0.0;
			      thompsonDepth = nEle * sigma_t * length;
			      if(thompsonDepth>=1.e-6)std::cerr<<"Thompson optical depth great > 1.e-6 \t"<<thompsonDepth<<"\n";
			      
			    }
			  jComp2 += dopp*jComp*(nPhotNu/phoNu)*(dNu*energyFactor);//integrate over nu
			  jComp = 0.0;
			}
		      //compton emission 
		      jComp2 *= (1./(c*h));
		      nPhotTheta += jComp2 * dTheta; //integration over photon theta
		      jComp2 = 0.0;
				      
		    }
		  nPhotPhi += nPhotTheta * dCosPhi;//integration over photon phi
		  nPhotTheta = 0.0;
		  
		}
	      iNu = mc2*c*phoEpsilonS*nPhotPhi *length;
	      intrInu = pho->get3DdistElement(itScatNu, itElePhi, itEleTheta, intrinsic);
	      iNu += intrInu;
	      pho->writeArr3Trans(itScatNu, itElePhi, itEleTheta, iNu, face);
	      nPhotPhi = 0.0;
	      iNu = 0.0;
	      	       
	    }
	}
    }
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Compton::comptonKernel(double &outPhotE, double &inPhotE_ERF, 
			      double &yParam, double &eleGamma)
{
  double kernel = yParam + (1./yParam) - ((2.*outPhotE)/(eleGamma*inPhotE_ERF*yParam))
    + sqr((outPhotE)/(eleGamma*inPhotE_ERF*yParam));
  return kernel;
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
// double Compton::comptoncsHeadOn(double &inPhotE_ERF, double &yParam, 
// 				double &eleGamma)
// {
//   double lowerLimit = 1./(1.+(2.*inPhotE_ERF));
//   double cs = 0.0;
  
//   cs = (physcon.pi*sqr(physcon.r_e)/(eleGamma*pow(inPhotE_ERF,3)))*
//     (yParam*sqr(inPhotE_ERF)+1.+(2.*inPhotE_ERF)+(1./yParam)*
//      (sqr(inPhotE_ERF)-(2.*inPhotE_ERF)-2.)+(1./sqr(yParam)));
//   return cs;
// }
