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

#ifndef CELL_HH
#define CELL_HH

#include <string>
#include <sstream>

#include "electrons.hpp"
#include "photons.hpp"
#include "hexa.hpp"
#include "synchrotron.hpp"
#include "compton.hpp"
#include "libs/numerical.hpp"
#include "libs/funcobj.hpp"

/*! Cell class for creating a single unit of volume*/
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
class Cell
{
  
public:
    
  Cell(const double &, const std::vector<int> &); 	//the default constructor
  Cell(const Cell &c);
  ~Cell();
  //create electron and photon distributions
  void electroPhotoDistributions();
  
  //set cell coords (string)
  void cellCoords();
  //the BField vector
  void setBField(double&, double &, double &);
  
   //set cell neighbours
  void setNeighbour(unsigned &, Cell *);
  
  //output neighbours
  void outputNeighbours();
  
  Cell * retNeighbour(int &);
  
  //Memeber function to carry out various physical processes
  void physicalProcesses();
 
  void passIntrPhToTrans(Photons *);
  
  void setPhotDirections();
  
  int faceIntersect(double &, double &, double &);
  
  int whichFace(double &, double &, double &);
  
  void passPhotonsToIntr();
  
  void createSynchArray();
    
  Photons * retPhot();
  Electrons * retElec();
  
  void passArrs(Photons *, int &);
  
  
  //return the doppler factor
  double dopplerFactor();
    
  //number of output faces for a cell
  int emptyFaceCount();
  
  Arr3D<int, double> * neighbourDist(int &);
  //set boundary cell i.e. not 6 neighbours
  
  double angleAberration(double);
   
  inline void setBoundaryCell()
  {
    boundary = true;
  }
  
  inline bool boundaryCell()
  {
    return boundary;
  }
  
  //return the cell position string
  inline std::string cellPos()
  {
    return position;
  }
  inline int cellX()
  {
    return coords[0];
  }
  inline int cellY()
  {
    return coords[1];
  }
  inline int cellZ()
  {
    return coords[2];
  }
  inline double B()
  {
    return BField;
  }
    
  inline double cellLength() const
  {
    return length;
  }
  
  inline double cellVolume() const
  { 
    return length * length * length;
  }
  
 
protected:
  
private:
  Electrons *ele;
  Photons *pho;
  Synchrotron *synch;
  Compton *compt;
  
  int grid;
  int agrid;
  
  double length;//cell size
  double BField;
  std::vector<double> BVectSpherical;  
  bool boundary;
  
  std::vector<int> coords;
  std::string position;
  Hexa<Cell *> neighbours;
  Numerical nums;
  GetIndex<double, double> getind;  
  Arr4D<int, double> *synchArray;
  
   
};
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....


#endif // CELL_HH
