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


#include "output.hpp"

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
Output::Output()
{
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
Output::~Output()
{
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
//writes electron file for every pitch angle
void Output::writeElectronFile(Electrons *ePtr, std::string &appen)
{
  std::vector<double>::const_iterator a, p;
  int a1 = 0, b1 =0;
  double tot1 = 0.0, eGamPhi = 0.0;
  std::string filename;
  std::stringstream ss;

  for (b1 = 0, p = ePtr->phiBegin(); b1 < ePtr->nSize();
       ++b1, ++p)
    {
      ss << std::setfill ('0')<< std::setw(3) << b1;
      filename = "elec_";
      filename += appen;
      filename += "_";
      filename += ss.str();
      filename += ".dat";
      std::ofstream elFile(filename.c_str());

      for (a = ePtr->xBegin(), a1=0; 
	   a != ePtr->xEnd(); ++a, ++a1)
	{
	  eGamPhi = ePtr->arrayElement2to1D(a1,b1);
	  elFile<<*a<<"\t"<<*p<<"\t"<<eGamPhi<<std::endl;
	}
      ss.str("");
    }
  
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
//writes photon file for every pitch angle
void Output::writePhotonFile(Photons *phPtr, std::string &appen)
{
  std::vector<double>::const_iterator a, p;
  int a1 = 0, b1 =0;
  double tot1 = 0.0, iNuPhi = 0.0;
  std::string filename;
  std::string phiFile = "phiFile.dat";
  std::ofstream phi(phiFile.c_str());
  //phi<<"number"<<"\t"<<"angle"<<std::endl;
  std::stringstream ss;
  
  for (b1 = 0, p = phPtr->phiBegin(); b1 < phPtr->nSize();
       ++b1, ++p)
    {
      //add zeros at the beginning of file name
      ss << std::setfill ('0')<< std::setw(3) << b1;
      filename = "phot_";
      filename += appen;
      filename += "_";
      filename += ss.str();
      filename += ".dat";
      std::ofstream nuFile(filename.c_str());

      for (a = phPtr->xBegin(), a1=0; 
	   a != phPtr->xEnd(); ++a, ++a1)
	{
	  iNuPhi = phPtr->arrayElement2to1D(a1,b1);
	  if (iNuPhi < 1.e-50)
	    {
	      iNuPhi = 0.0;
	    }
	  nuFile<<*a<<"\t"<<*p<<"\t"<<iNuPhi<<std::endl;
	    }
      phi<<b1<<"\t"<<*p<<std::endl;
      ss.str("");
      
    }
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
//writes photon and electron file for every boundary cell
void Output::writePhotElecOutput(std::vector<Cell *> *boundary)
{
  std::vector<Cell *>::const_iterator cell, cell2;
  int outputDist = 7, intrinsic  =6;
  int i,j,k;
  double tot1 = 0.0, iNu = 0.0, phi=0.0, theta=0.0, eDist=0.0;
  std::string filename, elFilename, cPos;
  Photons *ph;
  Electrons *el;
  std::stringstream ssNu, ssEl, ssNuBlock, ssElBlock;
  std::vector<double>::const_iterator nu, gam;
  
  std::cerr<<"Writing output files for boundary cells ";
  for (cell = boundary->begin(); cell != boundary->end();
       ++cell)
    {
      ph = (*cell)->retPhot();
      el = (*cell)->retElec();
      cPos = (*cell)->cellPos();
      //add zeros at the beginning of file name
      filename = "phot_cell_";
      elFilename = "el_cell_";
      filename += cPos;
      elFilename += cPos;
      filename += ".dat";
      elFilename += ".dat";
        
      std::ofstream nuFile(filename.c_str());
      std::ofstream elFile(elFilename.c_str());
      
      for (k = 0; k < ph->gridSize(); ++k)
	{
	  theta = ph->thetaElement(k);
	  for (j = 0; j < ph->agridSize(); ++j)
	    {
	      phi = ph->phiElement(j);
	      
	      //nuFile<<"theta="<<theta<<"\t"<<"phi="<<phi<<"\n";
	      //elFile<<"theta="<<theta<<"\t"<<"phi="<<phi<<"\n";
	      
	      for (nu = ph->xBegin(), gam = el->xBegin(), i=0; 
		   nu != ph->xEnd(); ++nu, ++gam, ++i)
		{
		  iNu = ph->get3DdistElement(i,j,k,outputDist);
		  eDist = el->get3DdistElement(i,j,k,intrinsic);
		  
		  // if (iNu < 1.e-50)
		  //   {
		  //     iNu = 0.0;
		  //   }
		  //using stringstream to store large blocks
		  // then output the whole block to a file; faster than simple << operator
		  //std::cout<<theta*(57.2957795)<<"\t"<<phi*(57.2957795)<<"\t"<<*nu<<"\t"<<iNu<<"\n";
		  ssNu<<theta*(57.2957795)<<"\t"<<phi*(57.2957795)<<"\t"<<*nu<<"\t"<<iNu<<"\n";
		  ssEl<<theta*(57.2957795)<<"\t"<<phi*(57.2957795)<<"\t"<<*gam<<"\t"<<eDist<<"\n";
		  //nuFile<<ssNu.rdbuf();
		  //elFile<<ssEl.rdbuf();
		  //nuFile.write(ssNu.str().c_str(), ssNu.str().length());
		  //elFile.write(ssEl.str().c_str(), ssNu.str().length());
		  //ssNu.str("");
		  //ssEl.str("");
		  
		}
	      ssNu<<"\n";
	      ssEl<<"\n";
	      //nuFile<<"\n";
	      //elFile<<"\n";
	      
	    }
	  ssNu<<"\n";
	  ssEl<<"\n";
	  nuFile.write(ssNu.str().c_str(), ssNu.str().length());
	  elFile.write(ssEl.str().c_str(), ssNu.str().length());
	  ssNu.str("");
	  ssEl.str("");
	  
	}
      ssNu.str("");
      ssEl.str("");
      theta=0.0;
      phi=0.0;
      iNu=0.0;
      
      nuFile.close();
      elFile.close();
      
   
            
    }
  // std::string nuF="nuFile.dat", gaF="gammaFile.dat", thetPhiF="thetPhi.dat";
  // cell2 = boundary->begin();
  // ph = (*cell2)->retPhot();
  // el = (*cell2)->retElec();
  // ph->outputNuRange(nuF);
  // el->outputGammaRange(gaF);
  // el->outputThetaPhi(thetPhiF);
  std::cerr<<"\n";
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Output::writeResArray(Arr3D<int, double> *resArr, 
			   std::vector<Cell *> *boundary, 
			   double &sourceDist, double &volArea, double &lFactor, 
			   double &redshift, int &tstep)
{
  std::string tSt = tToString(tstep);
  std::cerr<<"Writing results Array for all the boundary cells...";
  std::string filename = "output_"+tSt+".dat";
  std::string filename2 = "output_ave"+tSt+".dat";
  
  
  outputArray(filename, boundary, resArr, redshift, lFactor, volArea, sourceDist);
  outputArrayAverage(filename2, boundary, resArr, redshift, lFactor, volArea, sourceDist);
  std::cerr<<"...done"<<"\n";
      
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
//output array
void Output::outputArray(std::string &filename, 
			 std::vector<Cell *> *boundary, 
			 Arr3D<int, double> *resArr, double &redshift,
			 double &LF, double &VA, double &dist)
{
  std::vector<Cell *>::const_iterator cell = boundary->begin();
  std::stringstream ssNu;
  std::ofstream nuFile(filename.c_str());
 
  std::vector<double>::const_iterator nu;
  Photons *ph;
  ph = (*cell)->retPhot();
  int gridS = ph->gridSize();
  int agridS = ph->agridSize();
  double freqShift=0.0;
  double tot1 = 0.0, iNu = 0.0, phi=0.0, theta=0.0;
  int i, j, k;
  double val = 0.0, doppVal=0.0, dFlux=0.0;
  
  for (k = 0; k < agridS; ++k)
    {
      theta = ph->thetaElement(k);	  
          
      for (j = 0; j < agridS; ++j)
	{
	  phi = ph->phiElement(j);
	  doppVal = bulkDopplerFactor(phi, LF);
	  dFlux = dopplerFlux(doppVal, dist, VA);
	  freqShift = doppVal/(1.+redshift);
	  
	  for (i = 0, nu = ph->xBegin(); i < gridS; ++i, ++nu)	  
	    {
	      
	      iNu = resArr->getElement(i,j,k);
	      // ssNu<<((*cell)->angleAberration(theta))*(57.2957795)
	      // 	  <<"\t"<<((*cell)->angleAberration(phi))*(57.2957795)<<"\t"
	      // 	  <<*nu*freqShift<<"\t"<<iNu * dF<<"\n";
	      ssNu<<(theta)*(57.2957795)
		  <<"\t"<<((*cell)->angleAberration(phi))*(57.2957795)<<"\t"
		  <<*nu*freqShift<<"\t"<<iNu * dFlux<<"\n";
		
	    }
	  ssNu<<"\n\n";
	}
      ssNu<<"\n";
      nuFile.write(ssNu.str().c_str(), ssNu.str().length());
      ssNu.str("");
    }
  nuFile.close();

 
  
  std::cerr<<"\n";
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
//output array, but average over all the angles
//just a frequency vs. flux output from all the cells
void Output::outputArrayAverage(std::string &filename, 
				std::vector<Cell *> *boundary, 
				Arr3D<int, double> *resArr, double &redshift,
				double &LF, double &VA, double &dist)
{
  std::vector<Cell *>::const_iterator cell = boundary->begin();
  std::stringstream ssNu;
  std::ofstream nuFile(filename.c_str());
 
  std::vector<double>::const_iterator nu;
  Photons *ph;
  ph = (*cell)->retPhot();
  int gridS = ph->gridSize();
  int agridS = ph->agridSize();
  double freqShift = 0.0;
  double tot1 = 0.0, iNu = 0.0, phi=0.0, theta=0.0;
  int i, j, k;
  double val = 0.0, doppVal=0.0, dFlux=0.0;
 
  for (i = 0, nu = ph->xBegin(); i < gridS; ++i, ++nu)	  
    {
      for (k = 0; k < agridS; ++k)
	{
	  theta = ph->thetaElement(k);   

	  for (j = 0; j < agridS; ++j)
	    {
	      phi = ph->phiElement(j);
	      doppVal = bulkDopplerFactor(phi, LF);
	      dFlux = dopplerFlux(doppVal, dist, VA);
	      freqShift = doppVal/(1.+redshift);
	      iNu += resArr->getElement(i,j,k);
	    }
	}
      nuFile<<*nu*freqShift<<"\t"<<iNu * dFlux<<"\n";
      //nuFile.write(ssNu.str().c_str(), ssNu.str().length());
      iNu = 0.0;
      
    }
  nuFile.close();

 
  
  std::cerr<<"\n";
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Output::bulkDopplerFactor(double &source_angle, double &LF)
{
  //in comoving frame
  double lorentzF = LF;
  double beta = sqrt(1. - (1./sqr(lorentzF)));
  double viewTheta = source_angle;//radians
  double dopplerF = (lorentzF * (1+beta*(cos(viewTheta))));
  return dopplerF;
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Output::dopplerFlux(double &df, double &dist, double &volArea)
{
  double sourceDist = dist;
  double denom = 4. * physcon.pi * sqr(sourceDist);
  double val = (volArea/denom)*pow(df,3);
  return val;
      
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

