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

#ifndef HEXA_HH
#define HEXA_HH

//!Template class for defining directions in a cubic cell

#include <vector>
//Template class for creating a class with 6 directions
//It can be initiated with any type

template <typename T>
class Hexa : public std::vector<T> {
public:
  Hexa(const T &t) : std::vector<T>(6, t) 
  {}
  
  Hexa(const T &t1, const T &t2, const T &t3, const T &t4, const T &t5,
       const T &t6) : std::vector<T>() 
  { push_back(t1); push_back(t2); push_back(t3);
    push_back(t4); push_back(t5); push_back(t6); }
  
  //class destructor
  ~Hexa()
  {}
  
  //class member functions
  T &hexFront() { return *(this->begin()); }
  T &hexBack() { return *(this->begin()+1); }
  T &hexTop() { return *(this->begin()+2); }
  T &hexBottom() { return *(this->begin()+3); }
  T &hexRight() { return *(this->begin()+4); }
  T &hexLeft() { return *(this->begin()+5); }

}; 

//class for face vertex coordinates
template <typename T>
class ThreePs
{
public:
  ThreePs(const T &t1, const T &t2, const T &t3) :
    x(t1), y(t2), z(t3)
  {}
  
  T &xP() { return x; }
  T &yP() { return y; }
  T &zP() { return z; }

private:
  T x;
  T y;
  T z;
  
};

//class for storing 4 vertices of a square
template <typename T>
class Quaternary : std::vector<T>
{
public:
  Quaternary(const T &t1, const T &t2, const T &t3,
	  const T &t4): std::vector<T>()
  {
    push_back(t1); push_back(t2); push_back(t3);
    push_back(t4); 
  }
  ~Quaternary()
  {}
  
  //vertex 1, bottom left when looking at x-y plane
  //then counted clockwise
  T &quat1() { return *(this->begin()); }
  T &quat2() { return *(this->begin()+1); }
  T &quat3() { return *(this->begin()+2); }
  T &quat4() { return *(this->begin()+3); }

};

#endif // HEXA_HH
