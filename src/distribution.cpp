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
#include <cmath>
#include <fstream>

#include "distribution.hpp"
#include "libs/physcon.hpp"

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
/*! Distribution class for creating a generic container (std::vector<double>
  for the lepton and photon distributions*/
Distribution::Distribution(const int &size, const int &angleSize, 
			   const double &min, const double &max): grid(size),
								  angGrid(angleSize),
								  lower(min), 
								  upper(max)
					       
{
  //creating the log grid for lepton and photon distns.
  createLogGrid(grid, angGrid, lower, upper);
  
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
//copy constructor
// Distribution::Distribution(const Distribution *p)
//   :xPtr(p->xPtr), dxPtr(p->dxPtr),  nPtr(p->nPtr)   
// {}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
//!Class destructor deleting the containers
Distribution::~Distribution()
{
  // delete the vectors 
  delete xPtr, dxPtr, nPtr;/* phi, nPtrPhi;*/
  delete phi, theta, arrCube;
  delete trans0, trans1, trans2, trans3,
    trans4, trans5, outputDist, direction;  
  // Delete the 2D array for angle dependent distributions
  delete [] nPtrPhi;
  
 }
 //.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
 //!Creates the log grid as well 2D array for angle dependent distributions
 void Distribution::createLogGrid(const int &gridSize,
				  const int &agridSize,
				 const double &lowerBound, 
				 const double &upperBound)
 {
   //The number of "particles" vector
   nPtr = new std::vector<double>(gridSize,0.0);
   //x-grid
   xPtr = new std::vector<double>;
   //dx-grid
   dxPtr = new std::vector<double>;
   //Vector of a vector for phi (the angle between B and e^- velcity)
   phi = new std::vector<double>;
   theta = new std::vector<double>;
   arrCube = new Arr3D<int, double>(gridSize, agridSize, agridSize, 0.0);
   outputDist = new Arr3D<int, double>(gridSize, agridSize, agridSize, 0.0);
   
   direction = new Arr3D<int, int>(gridSize, agridSize, agridSize, 0);
   trans0 = new Arr3D<int, double>(gridSize, agridSize, agridSize, 0.0);
   trans1 = new Arr3D<int, double>(gridSize, agridSize, agridSize, 0.0);
   trans2 = new Arr3D<int, double>(gridSize, agridSize, agridSize, 0.0);
   trans3 = new Arr3D<int, double>(gridSize, agridSize, agridSize, 0.0);
   trans4 = new Arr3D<int, double>(gridSize, agridSize, agridSize, 0.0);
   trans5 = new Arr3D<int, double>(gridSize, agridSize, agridSize, 0.0);
   //array[rows][cols]
   //access array[i][j] by using array[j+i*Max_j]
   nPtrPhi = new double[gridSize * gridSize]; 
    
   //create grid for gamma, nu, and phi 
   double abscissa = lowerBound;
   double r = pow(upperBound / lowerBound, 1.0 / gridSize);
   double r2 = sqrt(r);
   double al = 0.01, al2 = 0.01;//-1. * (physcon.pi); //avoid singularity sin(0) = 0
  
   r -= 1.0;
     
   for (int i = 1; i <= gridSize; ++i) 
     {
       dxPtr->push_back(r * abscissa);
       xPtr->push_back(r2 * abscissa); // Calc log-space mid-point for grid
       abscissa += (r * abscissa);  // Add grid width to abscissa to find
      }
   for (int j = 1; j <= agridSize; ++j)
     {
       phi->push_back(al); //phi angle
       theta->push_back(al2); //azimuthal angle
       al += (physcon.pi-0.01)/(agridSize-1);
       al2 += (2.*physcon.pi-0.01)/(agridSize-1);
       //if changing the range, check on function Cell::newPhotCoords
     }
   
 }

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Distribution::differential(std::vector<double>::const_iterator pos) const
{
  double diffVal = 0.0;
  
  if (*pos == *nPtr->begin())
    {
      diffVal = (*(pos+2) - *pos)/(2.0*(*dxPtr)[1]);
    }
  else if (*pos == (*nPtr)[(nPtr->size())-1])
    {
      diffVal = (*pos - *(pos-2))/(2.0*(*dxPtr)[(nPtr->size()-2)]);
    }
  else if (*pos != *nPtr->begin() && *pos != (*nPtr)[(nPtr->size())-1])
    {
      diffVal = (*(pos+1) - *(pos-1))/(2.0*(*dxPtr)[pos-(nPtr->begin())]);
    }
  return diffVal;
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Distribution::differentialPhi(int &i, int &j) const
{
  double diffVal;
  if (i == 0)
    {
      diffVal = ((nPtrPhi[j+((i+2)*grid)]) - nPtrPhi[j+(i*grid)])/(2.0*(*dxPtr)[1]);
    }
  else if (i == (nPtr->size())-1)
    {
      diffVal = (nPtrPhi[j+(i*grid)] - nPtrPhi[j+((i-2)*grid)])/(2.0*(*dxPtr)[i-1]);
    }
  else 
    {
      diffVal = (nPtrPhi[j+((i+1)*grid)] - nPtrPhi[j+((i-1)*grid)])/(2.0*(*dxPtr)[i]);
    }
  return diffVal;
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Distribution::differentialPhiTheta(int &i, int &j, int &k, int &distId) const
{
  double diffVal;
  int i2 = i, i3 = i;
  
  if (i == 0)
    {
      i2 = i+2;
      diffVal = ((this->get3DdistElement(i2, j, k, distId)) - (this->get3DdistElement(i, j, k, distId)))/
	(2.0 *(*dxPtr)[1]);
    }
  else if (i == (nPtr->size())-1)
    {
      i2 = i-2;
      diffVal = ((this->get3DdistElement(i, j, k, distId)) - (this->get3DdistElement(i2, j, k, distId)))/
	(2.0 *(*dxPtr)[i-1]);
    }
  else 
    {
      i2 = i+1;
      i3 = i-1;
      diffVal = ((this->get3DdistElement(i2, j, k, distId)) - (this->get3DdistElement(i3, j, k, distId)))/
	(2.0 *(*dxPtr)[i]);
    }
  return diffVal;
}


//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Distribution::writeArr3Trans(int &i, int &j, int &k, double &val, 
				 int &face)
{
  if(face == 0)
    {
      trans0->putElement(i,j,k,val);
    }
  else if(face == 1)
    {
      trans1->putElement(i,j,k,val);
    }
  else if(face == 2)
    {
      trans2->putElement(i,j,k,val);
    }
  else if(face == 3)
    {
      trans3->putElement(i,j,k,val);
    }
  else if(face == 4)
    {
      trans4->putElement(i,j,k,val);
    }
  else if(face == 5)
    {
      trans5->putElement(i,j,k,val);
    }
   else if(face == 6)
    {
      arrCube->putElement(i,j,k,val);
    }
   else if(face == 7)
    {
      outputDist->putElement(i,j,k,val);
    }
  
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Distribution::get3DdistElement(int &i, int &j, int &k, int &id) const
{
  if(id == 0)
      {
	return (this->trans0)->getElement(i,j,k);
      }
  else if(id == 1)
      {
	return (this->trans1)->getElement(i,j,k);
      }
  else if(id == 2)
      {
	return (this->trans2)->getElement(i,j,k);
      }
  else if(id == 3)
      {
	return (this->trans3)->getElement(i,j,k);
      }
  else if(id == 4)
      {
	return (this->trans4)->getElement(i,j,k);
      }
  else if(id == 5)
      {
	return (this->trans5)->getElement(i,j,k);
      }
  else if(id == 6)
      {
	return (this->arrCube)->getElement(i,j,k);
      }
  else if(id == 7)
      {
	return (this->outputDist)->getElement(i,j,k);
      }
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
Arr3D<int, double> * Distribution::transiting(int &face)
{
  if(face == 0)
    {
      return this->trans0;
    }
  else if(face == 1)
    {      
      return this->trans1;
    }
  else if(face == 2)
    {
      return this->trans2;
    }
  else if(face == 3)
    {
      return this->trans3;
    }
  else if(face == 4)
    {
      return this->trans4;
    }
  else if(face == 5)
    {
      return this->trans5;
    }
  else if(face == 6)//returns the intrinsic photon distribution
    {
      return this->arrCube;
    }
  else if(face == 7)
    {
      return this->outputDist;
    }
  
}
  
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Distribution::outputThetaPhi(std::string &filename)
{
  std::vector<double>::const_iterator al, th;
  std::ofstream phiThetaFile(filename.c_str());
  
  for (al = phi->begin(), th = theta->begin()
	 ; al != phi->end(); ++al, ++th)
    {
      phiThetaFile<<*th*(57.2957795)<<"\t"<<*al*(57.2957795)<<"\n";
    }
  phiThetaFile.close();
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Distribution::resetArray(int &face)
{
  if(face == 0)
    {
      trans0->resetDoubArray();
    }
  else if(face == 1)
    {
      trans1->resetDoubArray();
    }
  else if(face == 2)
    {
      trans2->resetDoubArray();
    }
  else if(face == 3)
    {
      trans3->resetDoubArray();
    }
  else if(face == 4)
    {
      trans4->resetDoubArray();
    }
  else if(face == 5)
    {
      trans5->resetDoubArray();
    }
  else if(face == 6)
    {
      arrCube->resetDoubArray();
    }
  else if(face == 7)
    {
      outputDist->resetDoubArray();
    }
  else if(face == 8)
    {
      direction->resetIntArray();
    }
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Distribution::DistTotal(int distID)
{
 
  
  //get the relevant Array
  Arr3D<int, double> * arr = this->transiting(distID);
  int i, size = arr->arrSize();
  double element = 0.0, total = 0.0;
  
  for (i = 0; i < size-1; ++i)
    {
      element = arr->getElementSingle(i);
      total += element;
    }
  return total;
    
}



// int Distribution::getIndex(std::vector<double> *container, double &val)
// {
//   std::vector<double>::iterator lower;
  
//   lower = std::lower_bound(container->begin(), container->end(),
// 			   val);
//   if (lower == container->end())
//     {
//       return container->size();  
//     }
//   else
//     {
//       return int(lower - container->begin());
//     }
        
// }

// int Distribution::getIndex(std::vector<double> *container, double &val)
// {
//   std::vector<double>::iterator findit;
  
//   findit = std::find_if(container->begin(), container->end(),
// 			FindWithDouble(val));
//   if (findit == container->end())
//     {
//       return container->size();  
//     }
//   else
//     {
//       return int(findit - container->begin());
//     }
// }
