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

#ifndef DISTRIBUTION_HH
#define DISTRIBUTION_HH

#include <vector>
#include <iterator>
#include <algorithm>
#include <iostream>
#include "libs/funcobj.hpp"

//! Distributions class
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
class Distribution
{
  
public:
  //the default constructor
  Distribution(const int&, const int&, const double &, const double &);
  //Distribution(const Distribution & c);
  //Distribution(const Distribution *c);
  ~Distribution();
  
  // Return an iterator for modifying the grids
  std::vector<double>::iterator modBegin()
  { 
    return nPtr->begin();
  }  
  std::vector<double>::iterator modEnd()
  { 
    return nPtr->end(); 
  }
  //iterator to return number of lepton or photons
  std::vector<double>::const_iterator nBegin() const
  { 
    return nPtr->begin(); 
  }
  std::vector<double>::const_iterator nEnd() const
  { 
    return nPtr->end(); 
  }
  //iterator to return the grid values
  std::vector<double>::const_iterator xBegin() const
  { 
    return xPtr->begin(); 
  }
  std::vector<double>::const_iterator xEnd() const
  { 
    return xPtr->end(); 
  }

  double xVal(int &i)
  {
    return (*xPtr)[i];
  }
  
  //iterator to retrun dx for a given grid
  std::vector<double>::const_iterator dxBegin() const
  { 
    return dxPtr->begin(); 
  }     
  std::vector<double>::const_iterator dxEnd() const
  { 
    return dxPtr->end(); 
  }

  double dxVal(int &i)
  {
    return (*dxPtr)[i];
  }
  //return the nPtr itself
  std::vector<double> * nPReturn()
  {
    return nPtr;
  }
 
  //return the pitch angle iterator
  std::vector<double>::const_iterator phiBegin() const
  {
    return phi->begin();
  }
  std::vector<double>::const_iterator phiEnd() const
  {
    return phi->end();
  }

  std::vector<double>::const_iterator thetaBegin() const
  {
    return theta->begin();
  }
  std::vector<double>::const_iterator thetaEnd() const
  {
    return theta->end();
  }
  
  double phiElement(int &i)
  {
    return (*phi)[i];
  }
  double thetaElement(int &i)
  {
    return (*theta)[i];
  }
  
  std::vector<double> * retPhi() const
  {
    return phi;
  }

  std::vector<double> * retTheta() const
  {
    return theta;
  }
  
  //give the size of x-grid
  int nSize()
  {
    return xPtr->size();
  }
  //write to the 2D array of particles
  void write2to1DArray(int &i, int &j, double &val)
  {
    nPtrPhi[j+(i*grid)] = val;
  }

  //write to the 2D array of particles
  void write1DArray(int &i, double &val)
  {
    nPtrPhi[i] = val;
  }
  //return an element of the 1D array representing 2D array of particles
  double arrayElement2to1D(int &i, int &j)
  {
    return nPtrPhi[j+(i*grid)];
  }
  
  double arrayElement1D(int &i)
  {
    return nPtrPhi[i];
  }
  
  int nPtrPhiSize()
  {
    return grid * grid;
  }
  
  int gridSize()
  {
    return grid;
  }
  //the size of phi, theta grid
  int agridSize()
  {
    return angGrid;
  }
  
  //Array /////
  int arrSizeDist()
  {
    return arrCube->arrSize();
  }
     
  int getDirection(int &i, int &j, int &k)
  {
    return direction->getElement(i,j,k);
  }
  void writeDirection(int &i, int &j, int &k, int &val)
  {
    direction->putElement(i,j,k, val);
  }
     
  void addArray(Arr3D<int, double> *first, 
		Arr3D<int, double> *second)
  {
    first->addArray(second);
  }
  
  double get3DdistElement(int &, int &, int &, int &) const;
  //return the gradient at the given point on the grid/array
  double differential(std::vector<double>::const_iterator ) const;
  double differentialPhi(int &, int &) const; 
  double differentialPhiTheta(int &, int &, int &, int &) const; 
  //write transiting distributions
  void writeArr3Trans(int &, int &, int &, double &, 
		     int &);
 
  
  Arr3D<int, double> * transiting(int &);
  
  void outputThetaPhi(std::string &);
  void resetArray(int &);
  
  double DistTotal(int distID);
  
//int getIndex(std::vector<double> *, double &);
  

protected:
  void createLogGrid(const int &, const int &, const double &,
				  const double &);
  
  std::vector<double> *xPtr;
  std::vector<double> *dxPtr;
  std::vector<double> *nPtr;
  std::vector<double> *phi;
  std::vector<double> *theta;
  
  double *nPtrPhi;
  //3D array for storing distns as a function of phi and theta
  //the one 'intrinsic' distribution
  Arr3D<int, double> *arrCube, *outputDist;

  //array storing directions as funtion of angles and frequency/energy
  Arr3D<int, int> *direction;
  //Six 3D arrays for each face. These combine to form one intrinsic distrn
  Arr3D<int, double> *trans0;
  Arr3D<int, double> *trans1;
  Arr3D<int, double> *trans2;
  Arr3D<int, double> *trans3;
  Arr3D<int, double> *trans4;
  Arr3D<int, double> *trans5;
  
  //std::vector< std::vector<double> > nPtrPhi;
  
  
private:
  double lower;
  double upper;
  int grid;
  int angGrid;
  double zero;
    
};
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

#endif // DISTRIBUTION_HH
