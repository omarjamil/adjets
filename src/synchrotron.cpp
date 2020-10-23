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

#include <vector>
#include <iostream>
#include <cmath>

#include "libs/numerical.hpp"
#include "synchrotron.hpp"
#include "libs/funcobj.hpp"

#define sqr(a) ((a)*(a))

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

Synchrotron::Synchrotron()
{}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
Synchrotron::~Synchrotron()
{}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Synchrotron::synchPhiTheta(Electrons *ele, Photons *pho,
				const double &B, std::vector<double> &BVect, 
				const double &dl, const Arr4D<int, double> *synchArray)
{
  //std::cout<<"Synhcrotron Radiation...";
  //mostly Longair vol 2. treatment followed (SI units)
  Arr4D<int, double> synArray(*synchArray);
  
  double con1 = (sqrt(3.0) * pow(physcon.e,3)*B)
    /(4.*physcon.pi*physcon.m_e * physcon.c * physcon.epsilon_0);
  double con2 = (3*physcon.e*B)/(4.*physcon.pi*physcon.m_e);//nuCrit in SI
  //double con3 = (-1.*sqr(physcon.c))/(8.*physcon.pi);
  double con3 = (-1.)/(8.*physcon.pi*physcon.m_e);
  double gammaMin = *(ele->xBegin());
  double x, Pnu, jNu,  alpha, phi, thetaA, dndg=0.0, dnVar=0.0, ngPhi, nuSSA,
    absCoef=0.0, tauNu = 0.0, iNu=0.0, gamma, photNu, dgamma, intrInu, nuCritMin, sineAlpha;
  int n, ga, itNu, itAl, gIt=0,
    theta, face;
  int gridS = pho->gridSize();
  int agridS = pho->agridSize();
  int eDistId = 6;
  int phDistId, cellIntrinsic=6;
  for (itNu = 0; itNu <= gridS-1; ++itNu)
    {
      photNu = pho->xVal(itNu);// photon frequency
      for (itAl = 0; itAl <= agridS-1; ++itAl)
	{
	  phi = ele->phiElement(itAl);
	  
	  for (theta = 0; theta <= agridS-1; ++theta)
	    {
	      thetaA = ele->thetaElement(theta);
	      sineAlpha = sin(numerical.pitchAngle(phi,BVect[1],thetaA,BVect[2]));
	      if(sineAlpha < 0.0) sineAlpha *= -1.;
	      for (gIt = 0; gIt < gridS-1; ++gIt)
		{
		  gamma = ele->xVal(gIt);
		  //x2 as fast compared to synRad(x)
		  Pnu = synArray.getElement(itNu, itAl, theta, gIt);
		  ngPhi = ele->get3DdistElement(gIt, itAl, theta, eDistId);
		  dgamma = ele->dxVal(gIt);//dx of electron dist
		  //nuCritMin =  con2 * sqr(gammaMin) * sin(alpha);
		  //emission coefficient
		  jNu += Pnu * ngPhi * dgamma;
		  //absorption coefficient
		  dndg = ele->differentialPhiTheta(gIt, itAl, theta, eDistId);//differential
		  dnVar = dndg - ((2.0 * ngPhi) / gamma);
		  absCoef += dnVar * Pnu * dgamma;
		  Pnu = 0.0;
		}
	      jNu *= con1 * sineAlpha;//sin(alpha);
	      //std::cout<<con1<<"\t"<<sin(alpha)<<"\n";
	      //tauNu = con3*(1./sqr(photNu))*absCoef * dl;
	      tauNu = con1*sineAlpha*con3*(1./sqr(photNu))*absCoef * dl;
	      //nuSSA = pow(tauNu*sqr(photNu), 0.5);
	      //std::cout<<tauNu<<"\t"<<nuSSA<<"\n"
	  	      
	      if (tauNu > 1.e-10) //optically thick regime
		{
		  iNu = (jNu*dl/tauNu) * (1. - exp(-tauNu));
		}
	      else //optically thin
		{
		  iNu = jNu*dl;
		}
	      
	      // if(nuSSA > nuCritMin && photNu < nuCritMin)
	      // 	{
	      // 	  iNu = 0.0;
	      // 	}
	      
	      if(iNu <= 1.e-40) //arbitrary to eliminate really small values
		{
		  iNu = 0.0;
		}
	      //cell intrinsic photon distribution
	      intrInu = pho->get3DdistElement(itNu, itAl, theta, cellIntrinsic);  
	      intrInu *= exp(-tauNu);//absorption of the intrinsic
	      iNu += intrInu;
	    
	      face = 6;//pass to intrinsic and then to neighbours later
	      //(pho->transiting(face))->putElement(nuCount, itAl, theta, iNu);              	   
	      pho->writeArr3Trans(itNu, itAl, theta, iNu, face);
	      
	      
	      //pho->write1DArray(itNu, iNu);//write to photon distn.
	      jNu = 0.0, iNu=0.0, absCoef=0.0, tauNu = 0.0, nuSSA=0.0;
	    }
	}
    }
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Synchrotron::synRad(double &x)
/*   x :    energy normalized to the critical energy
      returns function value SynRadC   photon spectrum dn/dx
      (integral of modified 5/3 order Bessel function)
      principal: Chebyshev series see H.H.Umstaetter CERN/PS/SM/81-13 10-3-1981
      Also see BDSIM
*/

{ double synrad=0.;
  if(x>0. && x<800.) // otherwise result synrad remains 0
   { if(x<6.)
      { double a,b,z;
        z=x*x/16.-2.;
        b=          .00000000000000000012;
        a=z*b  +    .00000000000000000460;
        b=z*a-b+    .00000000000000031738;
        a=z*b-a+    .00000000000002004426;
        b=z*a-b+    .00000000000111455474;
        a=z*b-a+    .00000000005407460944;
        b=z*a-b+    .00000000226722011790;
        a=z*b-a+    .00000008125130371644;
        b=z*a-b+    .00000245751373955212;
        a=z*b-a+    .00006181256113829740;
        b=z*a-b+    .00127066381953661690;
        a=z*b-a+    .02091216799114667278;
        b=z*a-b+    .26880346058164526514;
        a=z*b-a+   2.61902183794862213818;
        b=z*a-b+  18.65250896865416256398;
        a=z*b-a+  92.95232665922707542088;
        b=z*a-b+ 308.15919413131586030542;
        a=z*b-a+ 644.86979658236221700714;
        double p;
        p=.5*z*a-b+  414.56543648832546975110;
        a=          .00000000000000000004;
        b=z*a+      .00000000000000000289;
        a=z*b-a+    .00000000000000019786;
        b=z*a-b+    .00000000000001196168;
        a=z*b-a+    .00000000000063427729;
        b=z*a-b+    .00000000002923635681;
        a=z*b-a+    .00000000115951672806;
        b=z*a-b+    .00000003910314748244;
        a=z*b-a+    .00000110599584794379;
        b=z*a-b+    .00002581451439721298;
        a=z*b-a+    .00048768692916240683;
        b=z*a-b+    .00728456195503504923;
        a=z*b-a+    .08357935463720537773;
        b=z*a-b+    .71031361199218887514;
        a=z*b-a+   4.26780261265492264837;
        b=z*a-b+  17.05540785795221885751;
        a=z*b-a+  41.83903486779678800040;
        double q;
        q=.5*z*a-b+28.41787374362784178164;
        double y;
        double const twothird=2./3.;
        y=pow(x,twothird);
        synrad=(p/y-q*y-1.)*1.81379936423421784215530788143;
      }
      else // 6 < x < 174
      { double a,b,z;
        z=20./x-2.;
        a=      .00000000000000000001;
        b=z*a  -.00000000000000000002;
        a=z*b-a+.00000000000000000006;
        b=z*a-b-.00000000000000000020;
        a=z*b-a+.00000000000000000066;
        b=z*a-b-.00000000000000000216;
        a=z*b-a+.00000000000000000721;
        b=z*a-b-.00000000000000002443;
        a=z*b-a+.00000000000000008441;
        b=z*a-b-.00000000000000029752;
        a=z*b-a+.00000000000000107116;
        b=z*a-b-.00000000000000394564;
        a=z*b-a+.00000000000001489474;
        b=z*a-b-.00000000000005773537;
        a=z*b-a+.00000000000023030657;
        b=z*a-b-.00000000000094784973;
        a=z*b-a+.00000000000403683207;
        b=z*a-b-.00000000001785432348;
        a=z*b-a+.00000000008235329314;
        b=z*a-b-.00000000039817923621;
        a=z*b-a+.00000000203088939238;
        b=z*a-b-.00000001101482369622;
        a=z*b-a+.00000006418902302372;
        b=z*a-b-.00000040756144386809;
        a=z*b-a+.00000287536465397527;
        b=z*a-b-.00002321251614543524;
        a=z*b-a+.00022505317277986004;
        b=z*a-b-.00287636803664026799;
        a=z*b-a+.06239591359332750793;
        double p;
        p=.5*z*a-b    +1.06552390798340693166;
        double const pihalf=physcon.pi/2.;
        synrad=p*sqrt(pihalf/x)/exp(x);
      }
   }
      
   return synrad;
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Synchrotron::synRadAsymptote(double &x)
{
  double synrad;
  if(x < 0.00001)
    {
      synrad = 4.*physcon.pi*(1./sqrt(3.))*pow((x/2.),1./3.)
	*numerical.gammalnFunction(1./3.);
    }
  else if(x > 1000.0)
    {
      synrad = sqrt(physcon.pi/2.)*exp(-1.*x)*sqrt(x);
    }

  return synrad;
    
}


