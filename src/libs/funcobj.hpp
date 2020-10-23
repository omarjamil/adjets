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

#ifndef FUNCOBJ_HH
#define FUNCOBJ_HH


//collection of function objects, templates etc.
// The following function objects are used as predicates for
// various finding/sorting/inserting algorithms provided by STL

#include <functional>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <iostream>
#include <algorithm>
#include <iterator>

#define sqr(a) ((a)*(a))
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
//3d array stored as 1D array
//array[row][col][..] = array[i][j][k]
//array[i][j][k] = array[j + i*max_j + k*max_i*max_j]
template <typename T, typename U>
class Arr3D 
{
public:
  //constructor
  Arr3D(const T &t, const U &u)
  {
    iS = t;
    jS = t;
    kS = t;
    arr = new U[(t*t*t)];
    allSize = iS * jS * kS;
    onet = true;
    
  }
 
  
  Arr3D(const Arr3D *a);
  
  Arr3D(const T &t1, const T &t2, const T &t3, const U &u)
  {
    iSize = t1;
    jSize = t2;
    kSize = t3;
    arr = new U[t1*t2*t3];
    initializeArr(u);
    
  }

  //copy constructor
  Arr3D (const Arr3D<T, U> &a)
  {
    iSize = (a.reti()); 
    jSize = (a.retj());
    kSize = (a.retk());
    arr = new U[iSize * jSize * kSize];
    copyArr(arr, (a.retarr()));
  }
  //destructor
  ~Arr3D()
  {
    delete [] arr;
    arr = NULL;
    
  }
 
  U *arr;
  unsigned int allSize, iSize, jSize, kSize, iS, jS, kS;
  bool onet;
  //for intialization
  void initializeArr(const U &u)
  {
    for(int i=0; i<(iSize * jSize * kSize); ++i)
      {
	(arr)[i] = u;
      }
  }
  

  //for copy constructor
  void copyArr(double arr[], double arr2[])
  {
    for(int i=0; i<(iSize * jSize * kSize); ++i)
      {
	(arr)[i] = (arr2)[i];
      }
  }
  
  U * retarr() const
  {
    return arr;
  }
  T reti() const
  {
    return iSize;
  }
  T retj() const
  {
    return jSize;
  }
  T retk() const
  {
    return kSize;
  }
  
  U getElementSingle(const T &s)
  {
    return this->arr[s];
  }
  
  U getElement(const T &i, const T &j, const T &k)
  {
    return this->arr[j + (i*jSize) + (k*iSize*jSize)];
  }
  
  U getElement(const T &i, const T &j, const T &k) const
  {
    return this->arr[j + (i*jSize) + (k*iSize*jSize)];
  }

  void putElement(const T &i, const T &j, const T &k, const U &e)
  {
    this->arr[j + (i*jSize) + (k*iSize*jSize)] = e;
  }
  
  void addElement(const T &i, const T &j, const T &k, const U &e)
  {
    U value = this->arr[j + (i*jSize) + (k*iSize*jSize)];
    U sum = value + e;
    this->arr[j + (i*jSize) + (k*iSize*jSize)] = sum;
  }

 
  T arrSize()
  {
    if(onet)
      {	
	return (iS * jS * kS);
      }
    else
      {
	return (iSize * jSize * kSize);
      }
    
  }
  //functions below not generic
  
  void resetDoubArray()
  {
    std::fill(&arr[0], &arr[iSize*jSize*kSize], 0.0);
  }
  void resetIntArray()
  {
    std::fill(&arr[0], &arr[iSize*jSize*kSize], 0);
  }
  void addArray(Arr3D<int, double> *second)
  {
    int i;
    double sum;
    int x = this->reti();
    int y = this->retj();
    int z = this->retk();
    for(i = 0; i < (x*y*z); ++i)
      {
	sum = this->arr[i] + (second->getElementSingle(i));
	this->arr[i] = sum;
      }
  }
};


//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
//Templates to convert a given type to string

template<typename T>
std::string tToString(T i)
{
  std::stringstream ss;
  std::string s;
  ss << i;
  ss >> s;
  return s;
}

//take 3 numbers and return one combined string
template<typename T>
std::string tToString3(T i, T j, T k)
{
  std::stringstream ss, tt, uu;
  std::string s, t, u, full;
  ss << i;
  tt << j;
  uu << k;
  ss >> s;
  tt >> t;
  uu >> u;
  full = s;
  full += t;
  full += u;
  return full;
}
//Fill can take an integer and fill zeros in front
template<typename T>
std::string tToStringFill(T t, int i)
{
  std::stringstream ss;
  std::string s;
  ss << std::setfill('0')<<std::setw(i)<<t;
  ss >> s;
  return s;
}

template<typename T>
std::vector<int> stringToVecInt(T s)
{
  int x, y, z;
  std::vector<int> coord;
  std::stringstream xs,ys,zs;
  std::string xx = s.substr(0,1);
  xs<<xx;
  xs>>x;
  coord.push_back(x);
  
  std::string yy = s.substr(1,1);
  ys<<yy;
  ys>>y;
  coord.push_back(y);
  
  std::string zz = s.substr(2,1);
  zs<<zz;
  zs>>z;
  coord.push_back(z);
  
  return coord;
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
template <typename T, typename U>
class GetIndex
{
   
public:
  GetIndex(){};
  ~GetIndex(){};
  

  //Get the nearest match
  int result(std::vector<T> *container, U val)
  {
    typename std::vector<T>::const_iterator it = std::lower_bound(container->begin(), 
				     container->end(), val);
    if (it == container->end()) 
      { 
	return container->size() - 1; 
      }
    else
      {
	return int(it-container->begin());
      }
    
  }
  
};


//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
//used to find a value in a container set up for a double
template <typename T>
class FindWithVal
{
  T t;
  T retVal;   
public:
  FindWithVal(T &s) : t(s){}
   

  bool operator()(T &val)
  {
    retVal = t-val;
    return (retVal < 0.015 && retVal > -0.1);
  }
};




  
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
// Function object structure that inherits from STL binary
// functions. Return operator () defined to return a result.
// Only a bool required as this is used with STL 'stable_sort' as 
// 'comp' function. Structure used as only a simple thing with only
// public members required, no initialization. No arguments being
// passed. Class can be used with public decl only.

// class IntToString
// {
//   int input;
  
// public:
//   IntToString(int &in): input(in)
//   {}
  
//   std::string operator()(int input) 
//   {
//     std::string line;
//     std::stringstream ss;
//     ss << std::setfill('0') << std::setw(3) << input;
//     ss >> line;
//     return line;
//   }

// };

// template<typename T>
// class GetIndex
// {
//   std::vector<T> *container;
//   T val;
 
  
// public:
//   GetIndex(std::vector<T> *con, T &v) : 
//     container(con), val(v) 
//   {
   
//   }

//   int operator()(std::vector<T> *con, T &v)
//   {
   
//     typename std::vector<T>::iterator findit;   
//     findit = std::find_if(container->begin(), container->end(),
// 			  FindWithVal(val));
//     if (findit == container->end())
//       {
// 	return container->size();  
//       }
//     else
//       {
// 	return int(findit - container->begin());
//       }
//   }
  
  
// };


//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
template< typename T >
struct delete_ptr : public std::unary_function<bool,T>
{
  bool operator()(T*pT) 
    const 
  { 
    delete pT; 
    return true; 
  }
};
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
//4d array stored as 1D array
//array[row][col][..] = array[i][j][k][l]
//array[i][j][k][l] = array[j + i*max_j + k*max_i*max_j + l*max_i*max_j*max_k]

template <typename T, typename U>
class Arr4D 
{
public:
  //constructor
  Arr4D(const T &t, const U &u)
  {
    iS = t;
    jS = t;
    kS = t;
    lS = t;
    arr = new U[(t*t*t*t)];
    allSize = iS * jS * kS * lS;
    onet = true;
  }
  Arr4D(const Arr4D *a);
  
  Arr4D(const T &t1, const T &t2, const T &t3, const T &t4, const U &u)
  {
    iSize = t1;
    jSize = t2;
    kSize = t3;
    lSize = t4;
    arr = new U[t1*t2*t3*t4];
    initializeArr(u);
    
  }

  //copy constructor
  Arr4D (const Arr4D<T, U> &a)
  {
    iSize = (a.reti()); 
    jSize = (a.retj());
    kSize = (a.retk());
    lSize = (a.retl());
    arr = new U[iSize * jSize * kSize * lSize];
    copyArr(arr, (a.retarr()));
  }
  //destructor
  ~Arr4D()
  {
    delete [] arr;
    arr = NULL;
    
  }
 
  U *arr;
  unsigned int allSize, iSize, jSize, kSize, lSize, iS, jS, kS, lS;
  bool onet;
  //for intialization
  void initializeArr(const U &u)
  {
    for(int i=0; i<(iSize * jSize * kSize * lSize); ++i)
      {
	(arr)[i] = u;
      }
  }
  

  //for copy constructor
  void copyArr(double arr[], double arr2[])
  {
    for(int i=0; i<(iSize * jSize * kSize * lSize); ++i)
      {
	(arr)[i] = (arr2)[i];
      }
  }
  
  U * retarr() const
  {
    return arr;
  }
  T reti() const
  {
    return iSize;
  }
  T retj() const
  {
    return jSize;
  }
  T retk() const
  {
    return kSize;
  }

  T retl() const
  {
    return lSize;
  }
  
  U getElementSingle(const T &s)
  {
    return this->arr[s];
  }
  
  U getElement(const T &i, const T &j, const T &k, const T &l)
  {
    return this->arr[j + (i*jSize) + (k*iSize*jSize) + (l*iSize*jSize*kSize)];
  }
  
  U getElement(const T &i, const T &j, const T &k, const T &l) const
  {
    return this->arr[j + (i*jSize) + (k*iSize*jSize) + (l*iSize*jSize*kSize)];
  }

  void putElement(const T &i, const T &j, const T &k, const T &l, const U &e)
  {
    this->arr[j + (i*jSize) + (k*iSize*jSize) + (l*iSize*jSize*kSize)] = e;
  }
  
  void addElement(const T &i, const T &j, const T &k, const T &l, const U &e)
  {
    U value = this->arr[j + (i*jSize) + (k*iSize*jSize) + (l*iSize*jSize*kSize)];
    U sum = value + e;
    this->arr[j + (i*jSize) + (k*iSize*jSize) + (l*iSize*jSize*kSize)] = sum;
  }

 
  T arrSize()
  {
    if(onet)
      {	
	return (iS * jS * kS * lS);
      }
    else
      {
	return (iSize * jSize * kSize * lSize);
      }
    
  }
  //functions below not generic
  
  void resetDoubArray()
  {
    std::fill(&arr[0], &arr[iSize*jSize*kSize*lSize], 0.0);
  }
  void resetIntArray()
  {
    std::fill(&arr[0], &arr[iSize*jSize*kSize*lSize], 0);
  }
  
  void addArray(Arr4D<int, double> *second)
  {
    int i;
    double sum;
    int x = this->reti();
    int y = this->retj();
    int z = this->retk();
    int a = this->retl();
    
    for(i = 0; i < (x*y*z*a); ++i)
      {
	sum = this->arr[i] + (second->getElementSingle(i));
	this->arr[i] = sum;
      }
  }
};




#endif // FUNCOBJ_HH
