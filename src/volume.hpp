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

#ifndef VOLUME_HH
#define VOLUME_HH

#include <vector>
#include <map>

#include "electrons.hpp"
#include "photons.hpp"
#include "synchrotron.hpp"
#include "output.hpp"
#include "cell.hpp"

//! The physical volume is created here
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
class Volume
{
  
public:
  Volume(const double &); 	//the default constructor
  Volume(const Volume &c);
  ~Volume();
  //all the processes are carried out in a cell
 
  
  //Create cells and populate the cell container
  void createVolume();
  //link the created cells
  void linkCells(std::vector<Cell *> *);
  //initialize the cells - eventually read from a file for every cell
  void initializeCells();
  
  //iterate through the cell container and carry out the physical processes
  void evolveCells();
  void loadBFile(std::string &, double &, double &, double &);
  
  void evolveVolume();
    
  void combineBoundaryOutput();
  
  inline double volArea()
  {
    return volumeArea;
  }
  
  inline void setVolArea(double volA)
  {
    volumeArea = volA;
  }
  
  double dopplerFlux(double &);
  double bulkDopplerFactor(double &);
  
protected:
    
private:
  friend class Output;
  Output *outp;
  Cell *cells;
  
  double BField;
  std::vector<double> BFieldVec;
  double cellSize;
  unsigned dir;
  int faceCount;
  double volumeArea;
  
  std::vector<Cell *> *cellContainer;
  std::vector<Cell *> *boundaryCells;  
  //store the position of a cell in the container
  //the key is the cell position (string form)
  std::map<std::string, int> cellMap;
  Arr3D<int, double> *resultsArray;
  //give the luminosity distance
  double NedWright();   
 
};
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....


#endif // VOLUME_HH
