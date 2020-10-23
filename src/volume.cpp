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


#include <iostream>
#include <fstream>
#include <cmath>

#include "volume.hpp"
#include "paraccess.hpp"
#include "libs/funcobj.hpp"
#include "photons.hpp"

#define sqr(a) ((a)*(a))

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
Volume::Volume(const double &size) : cellSize(size), volumeArea(0.0) 
{
  createVolume();
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
Volume::~Volume()
{
  for (std::vector<Cell *>::iterator it = cellContainer->begin(); 
       it != cellContainer->end(); ++it)
    {
      delete *it;
    }
  delete cellContainer;
  delete boundaryCells;
  delete outp;
  
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Volume::createVolume()
{
  cellContainer = new std::vector<Cell *>;
  
  const double size = cellSize;
  std::vector<int> position(3,0);
  std::string cPos, filename="cellVis.dat";
  std::stringstream ss;
  
  const int x = ext_x_cells;
  const int y = ext_y_cells;
  const int z = ext_z_cells;
  int i,j,k,s, cCount=0;
  std::cout<<"Creating cells... ";
  std::ofstream cellFile(filename.c_str());
  //loops to create the desired number of cells in the x, y, z directions
  for (i=0;i<x;++i)
    {
      
      for (j=0;j<y;++j)
	{
	 
	  for (k=0;k<z;++k)
	    {
	      ++cCount;
	      position[0]=i;
	      position[1]=j;
	      position[2]=k;
	      cells = new Cell(size, position);
	      cPos = cells->cellPos();
	      ss<<i<<"\t"<<j<<"\t"<<k<<"\t"<<cPos<<"\n";
	      cellFile.write(ss.str().c_str(), ss.str().length());
	      ss.str("");
	      
	      cellContainer->push_back(cells);
	      s=cellContainer->size()-1;//the number of cell in the container
	      cellMap[cPos]=s;
	    }
	}
    }
  cellFile.write(ss.str().c_str(), ss.str().length());
  std::cout<<"Created: "<<cCount<<" cells"<<std::endl;
  initializeCells();//initialize the B-Field from a file for now
  linkCells(cellContainer);
  
  
  
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Volume::linkCells(std::vector<Cell *> *cContainer)
{
  
  int xx, yy, zz, xp, xn, yp, yn, zp, zn, contPos;
  int xMax = ext_x_cells -1, yMax = ext_y_cells -1, zMax = ext_z_cells -1;
  std::string cellP;
  
  std::vector<Cell *>::const_iterator cIt;
   
  std::map<std::string, int>::iterator mapIt;
    
  Cell *bCells, *tempCell;
  boundaryCells = new std::vector<Cell *>; 
  bool edgeCell;
  
  std::cout<<"Linking cells... ";
  
  for(cIt = cContainer->begin(); cIt != cContainer->end(); ++cIt)
    {
      xx = (*cIt)->cellX();
      yy = (*cIt)->cellY();
      zz = (*cIt)->cellZ();
   
      //+ve x neighbour: FORWARD
      xp = xx + 1;
      if(xp <= xMax)
	{
	  dir = FORWARD;
	  cellP = tToString3(xp, yy, zz);
	  mapIt = cellMap.find(cellP);
	  contPos = mapIt->second;
	  tempCell = (*cContainer)[contPos];
	  (*cIt)->setNeighbour(dir, tempCell);
	}
      else
	{
	  bCells = *cIt;
	  edgeCell = bCells->boundaryCell();
	  if(!edgeCell)
	    {
	      boundaryCells->push_back(bCells);
	      bCells->setBoundaryCell();
	    }
	}
      //+ve x neighbour: BACK
      xn = xx -1;
      if (xn >= 0)
	{
	  dir = BACKWARD;
	  cellP = tToString3(xn, yy, zz);
	  mapIt = cellMap.find(cellP);
	  contPos = mapIt->second;
	  tempCell = (*cContainer)[contPos];
	  (*cIt)->setNeighbour(dir, tempCell);
	}
       else
	{
	  bCells = *cIt;
	  edgeCell = bCells->boundaryCell();
	  if(!edgeCell)
	    {
	      boundaryCells->push_back(bCells);
	      bCells->setBoundaryCell();
	    }
	}
      //+ve y: UP
      yp = yy + 1;
      if (yp <= yMax)
	{
	  dir = UP;
	  cellP = tToString3(xx, yp, zz);
	  mapIt = cellMap.find(cellP);
	  contPos = mapIt->second;
	  tempCell = (*cContainer)[contPos];
	  (*cIt)->setNeighbour(dir, tempCell);
	}
      else
	{
	  bCells = *cIt;
	  edgeCell = bCells->boundaryCell();
	  if(!edgeCell)
	    {
	      boundaryCells->push_back(bCells);
	      bCells->setBoundaryCell();
	    }
	}
      //-ve y: DOWN
      yn = yy - 1;
      if (yn >= 0)
	{
	  dir = DOWN;
	  cellP = tToString3(xx, yn, zz);
	  mapIt = cellMap.find(cellP);
	  contPos = mapIt->second;
	  tempCell = (*cContainer)[contPos];
	  (*cIt)->setNeighbour(dir, tempCell);
	}
       else
	{
	  bCells = *cIt;
	  edgeCell = bCells->boundaryCell();
	  if(!edgeCell)
	    {
	      boundaryCells->push_back(bCells);
	      bCells->setBoundaryCell();
	    }
	}
      //+ve z: RIGHT
      zp = zz + 1;
      if (zp <= zMax)
	{
	  dir = RIGHT;
	  cellP = tToString3(xx, yy, zp);
	  mapIt = cellMap.find(cellP);
	  contPos = mapIt->second;
	  tempCell = (*cContainer)[contPos];
	  (*cIt)->setNeighbour(dir, tempCell);
	}
      else
	{
	  bCells = *cIt;
	  edgeCell = bCells->boundaryCell();
	  if(!edgeCell)
	    {
	      boundaryCells->push_back(bCells);
	      bCells->setBoundaryCell();
	    }
	}
      //-ve z: LEFT
      zn = zz - 1;
      if (zn >= 0)
	{
	  dir = LEFT;
	  cellP = tToString3(xx, yy, zn);
	  mapIt = cellMap.find(cellP);
	  contPos = mapIt->second;
	  tempCell = (*cContainer)[contPos];
	  (*cIt)->setNeighbour(dir, tempCell);
	}
       else
	{
	  bCells = *cIt;
	  edgeCell = bCells->boundaryCell();
	  if(!edgeCell)
	    {
	      boundaryCells->push_back(bCells);
	      bCells->setBoundaryCell();
	    }
	}
    }
  std::cout<<"... done"<<std::endl;
  
  
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Volume::initializeCells()
{
  std::vector<Cell *>::iterator cellIt;
  std::string cPos;
  double BField, BTheta, BPhi;
  //loop through the cells in the x, y, z directions
 
  for (cellIt = cellContainer->begin(); 
       cellIt != cellContainer->end(); ++cellIt)
    {
      cPos = (*cellIt)->cellPos();
      loadBFile(cPos, BField, BTheta, BPhi);
      (*cellIt)->setBField(BField, BTheta, BPhi);
      (*cellIt)->setPhotDirections();
      (*cellIt)->createSynchArray();
      BField = 0.0;
     
    }
	 
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Volume::evolveVolume()
{
  std::vector<Cell *>::const_iterator cell = boundaryCells->begin();
  int tSteps = ext_t_steps;
  double tStepSize = ((*cell)->cellLength())/3.e8;
  std::cerr<<"Each time step is: "<<tStepSize<<" s"<<"\n";
  double redshift = ext_source_z;
  double LorentzFactor = ext_lorentz_factor;
  double volAr = 0.0;
  double inclination_angle = ext_source_angle;
  double doppVal = this->bulkDopplerFactor(inclination_angle);
  
  double sourceDist = this->NedWright() * 3.08568025e22;
  for(int i = 1; i <= tSteps; ++i)
    {
      std::cerr<<"Time step: "<<i<<"\n";
      this->evolveCells();
      this->combineBoundaryOutput();
      volAr = this->volArea();
      outp->writeResArray(resultsArray, boundaryCells, sourceDist, volAr, 
			  LorentzFactor, redshift, i);
           
    }

  //std::vector<Cell *>::const_iterator cell = boundaryCells->begin();
  Photons *ph;
  Electrons *el;
  ph = (*cell)->retPhot();
  el = (*cell)->retElec();  
  std::string nuF="nuFile.dat", gaF="gammaFile.dat", thetPhiF="thetPhi.dat",
    thetPhi="cellThetaPhi.dat";

  //ph->outputNuRange(nuF, redshift, doppVal);
  el->outputGammaRange(gaF);
  el->outputThetaPhi(thetPhiF);
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Volume::evolveCells()
{
  std::vector<Cell *>::const_iterator cellIt, cellIt2;
  int count = 1;
  outp = new Output();
  
  cellIt = cellContainer->begin();
  while(cellIt != cellContainer->end())
    {
      std::cerr<<"\r"<<"Physical processes in cell: "<<count<<" ";
      (*cellIt)->physicalProcesses();
      ++cellIt;        
      ++count;
    }
   
  count = 1;
  cellIt2 = cellContainer->begin();
  while(cellIt2 != cellContainer->end())
    {
      std::cerr<<"\r"<<"Passing photons to neighbours: cell "<<count<<" ";
     
      (*cellIt2)->passPhotonsToIntr();
  
      
      ++cellIt2;
      ++count;
    }
  std::cerr<<"\n";
  //outp->writePhotElecOutput(boundaryCells);

  // std::vector<Cell *>::const_iterator cellTest;
  // cellTest = cellContainer->begin();
  // Photons *phTest = (*cellTest)->retPhot();
  // double t = 1.00675391;
  // std::vector<double> *vecTest = phTest->retPhi();
  // int t2 = phTest->getIndex(vecTest, t);
  // int t3 = t2-1;
  // std::cout<<t<<"\t"<<t2<<"\t"<<phTest->phiElement(t2)<<"\t"<<phTest->phiElement(t3)<<std::endl;
  
    
  
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Volume::loadBFile(std::string &cellPos, 
		       double &B, double &BTheta, double &BPhi)
{
  std::string filename = ext_B_file;
  double B_theta_tmp, B_phi_tmp;
  
  std::ifstream BFile(filename.c_str());
  if(!BFile.good())
    {
      std::cerr<<"B-field file failed to load"<<std::endl;
    }
  std::string label, dummy;
  
   while(!BFile.eof())
     {
       BFile >> label;
       
       if(label == cellPos)
	 {
	   BFile >> dummy >> B >> B_theta_tmp >> B_phi_tmp;
	 }
     }
   BTheta = B_theta_tmp * 0.017453292;
   BPhi = B_phi_tmp * 0.017453292;
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Volume::combineBoundaryOutput()
{
  std::vector<Cell *>::const_iterator cellIt;
  int gridSize = ext_dist_grid;
  int agridSize = ext_angle_dist_grid;
  
  resultsArray = new Arr3D<int, double>(gridSize, agridSize, agridSize, 0.0);
  Arr3D<int, double> *cellOutput;
  Photons *phots;
  
  int output = 7;
  faceCount = 0;
  
  for (cellIt = boundaryCells->begin(); cellIt != boundaryCells->end(); 
       ++cellIt)
    {
      phots = (*cellIt)->retPhot();
      faceCount += (*cellIt)->emptyFaceCount();
      cellOutput = phots->transiting(output);
      resultsArray->addArray(cellOutput);
      phots->resetArray(output);
      
    }
  
  //volumeArea = faceCount * sqr(cellSize);
  this->setVolArea(faceCount * sqr(cellSize));
    
    
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Volume::dopplerFlux(double &df)
{
  //std::vector<Cell *>::const_iterator cellIt;
  //cellIt = boundaryCells->begin();
  //double df = (*cellIt)->dopplerFactor();
  //double df = this->bulkDopplerFactor();
  //double sourceDist = ext_source_d * 3.08568025e19;
 
  
  double sourceDist = this->NedWright() * 3.08568025e22;
  double denom = 4. * physcon.pi * sqr(sourceDist);
  double val = ((this->volArea())/denom)*pow(df,3);
  return val;
      
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

double Volume::NedWright()
{
  // Reference: Wright, E. L., 2006, PASP, 118, 1711.
  double z = ext_source_z;
  int n=1000;	// number of points in integrals
  double c = 299792.458; // velocity of light in km/sec
  double Tyr = 977.8; // coefficent for converting 1/H into Gyr
  double h0 = 0.7;
  double H0 = h0*100.;	// Hubble constant
  double Omega_M = 0.270;
  double Omega_Lambda = 0.730;
  
  // Densities
  double WM = Omega_M;		// Omega(matter)
  double WV = Omega_Lambda;	// Omega(vacuum) or lambda
  double h = H0/100.;	// H0/100
  double WR = 4.165e-5/(h*h);	// Omega(radiation), includes 3 massless neutrino species, T0 = 2.72528
  double WK = 1-WM-WR-WV;	// Omega curvaturve = 1-Omega(total)
  
  double a = 1.0;	// the scale factor of the Universe
  double az = 1.0/(1.+1.0*z);
  
  double age = 0;
  double adot;
  for (int i = 0; i < n; i++) 
    {
      a = az*((double)i+0.5)/(double)n;
      adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a));
      age = age + 1/adot;
    }
  double zage = az*age/n;
  // correction for annihilations of particles not present now like e+/e-
  // added 13-Aug-03 based on T_vs_t.f
  double lpz = log((1.+1.0*z))/log(10.0);
  double dzage = 0.;
  if (lpz >  7.500) dzage = 0.002 * (lpz -  7.500);
  if (lpz >  8.000) dzage = 0.014 * (lpz -  8.000) +  0.001;
  if (lpz >  8.500) dzage = 0.040 * (lpz -  8.500) +  0.008;
  if (lpz >  9.000) dzage = 0.020 * (lpz -  9.000) +  0.028;
  if (lpz >  9.500) dzage = 0.019 * (lpz -  9.500) +  0.039;
  if (lpz > 10.000) dzage = 0.048;
  if (lpz > 10.775) dzage = 0.035 * (lpz - 10.775) +  0.048;
  if (lpz > 11.851) dzage = 0.069 * (lpz - 11.851) +  0.086;
  if (lpz > 12.258) dzage = 0.461 * (lpz - 12.258) +  0.114;
  if (lpz > 12.382) dzage = 0.024 * (lpz - 12.382) +  0.171;
  if (lpz > 13.055) dzage = 0.013 * (lpz - 13.055) +  0.188;
  if (lpz > 14.081) dzage = 0.013 * (lpz - 14.081) +  0.201;
  if (lpz > 15.107) dzage = 0.214;
  zage = zage*pow(10.0,dzage);
  double zage_Gyr = (Tyr/H0)*zage;
  
  double DTT = 0.0;
  double DCMR = 0.0;
  // do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
  for (int i = 0; i < n; i++) 
    {
      a = az+(1.-az)*((double)i+0.5)/(double)n;
      adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a));
      DTT = DTT + 1./adot;
      DCMR = DCMR + 1./(a*adot);
    }
  
  DTT = (1-az)*DTT/(double)n;
  DCMR = (1-az)*DCMR/(double)n;
  
  age = DTT+zage;
  double age_Gyr = age*(Tyr/H0);
	
  
  double DTT_Gyr = (Tyr/H0)*DTT;
	
  double DCMR_Gyr = (Tyr/H0)*DCMR;
  double DCMR_Mpc = (c/H0)*DCMR;
  
  double DA, DA_Gyr, DL;
  
  // tangential comoving distance
  double ratio = 1.00;
  double x;
  double y;
  x = sqrt(fabs(WK))*DCMR;
  if (x > 0.1) 
    {
      if (WK>0) ratio = 0.5*(exp(x)-exp(-x))/x;
      else ratio = sin(x)/x;
      y = ratio*DCMR;
    }
  else	
    {
      y = x*x;
      // statement below fixed 13-Aug-03 to correct sign error in expansion
      if (WK < 0) y = -y;
      ratio = 1. + y/6. + y*y/120.;
      y= ratio*DCMR;
    }
  double DCMT = y;
  
  DA = az*DCMT;
  double DA_Mpc = (c/H0)*DA;
  double kpc_DA = DA_Mpc/206.264806;
  DA_Gyr = (Tyr/H0)*DA;
  DL = DA/(az*az);
  double DL_Mpc = (c/H0)*DL;
  double DL_Gyr = (Tyr/H0)*DL;
  
  // printf("Distance parameters: age_Gyr %e zage_Gyr %e DTT_Gyr %e DA_Mpc %e kpc_DA %e DL_Mpc %e \n", *age_Gyr,*zage_Gyr,*DTT_Gyr,*DA_Mpc,*kpc_DA,*DL_Mpc);
	
  return DL_Mpc;
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Volume::bulkDopplerFactor(double &source_angle)
{
  double lorentzF = ext_lorentz_factor;
  double beta = sqrt(1. - (1./sqr(lorentzF)));
  double viewTheta = source_angle;//radians
  double dopplerF = 1./(lorentzF * (1-beta*(cos(viewTheta))));
  return dopplerF;
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

