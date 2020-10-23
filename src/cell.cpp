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

#include <fstream>
#include "cell.hpp"
#include "paraccess.hpp"


#define sqr(a) ((a)*(a))

/*! Cell class code for creating a single unit of volume where all the physical
processes are carried out.*/
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
Cell::Cell(const double &d1, const std::vector<int> &pos):
  length(d1), BField(0.0), boundary(false),
  coords(pos), neighbours(0), BVectSpherical(3,0.0)
{
  cellCoords();
  electroPhotoDistributions();
  synch = new Synchrotron();
  compt = new Compton();
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
Cell::~Cell()
{
  delete ele;
  delete pho;
  delete synch;
  delete synchArray;
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

void Cell::electroPhotoDistributions()
{
  const int phSize = ext_dist_grid;
  const double phMin = ext_ph_nu_min;
  const double phMax = ext_ph_nu_max;
  const int eleSize = ext_dist_grid;
  const double eleMin = ext_e_gamma_min;
  const double eleMax =  ext_e_gamma_max;
  const double eleBrk = ext_e_gamma_brk;
  const double kT = ext_kT_e;
  const double neTh = ext_therm_norm;
  const double neP = ext_pLaw_norm;
  const double p = ext_pLaw_index;
  const double p1 = ext_pLaw_index1;
  const double r1 = ext_cell_size;
  const double volume = cellVolume();
  const double eleEneDens = ext_elec_ke_dens;
  const bool use_norm = ext_use_norm_bool;
  const int asize = ext_angle_dist_grid;
  
  ele =  new Electrons(eleSize, asize, eleMin, eleMax);
  ele->thermal(neTh, kT);
  if(use_norm)
    {
      ele->powerLaw(neP, p, eleMin, volume);
    }
  else
    {
      ele->powerLawEnDen(eleEneDens, p, p1, eleMin, eleMax, eleBrk, volume);
    }
  
  //elePtr->powerLawPhi();
  ele->powerLawPhiTheta();
  
  pho = new Photons(phSize, asize, phMin, phMax);
}


//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Cell::cellCoords()
{
  std::string x = tToString(coords[0]);
  position = x;
  std::string y = tToString(coords[1]);
  position += y;
  std::string z = tToString(coords[2]);
  position += z;
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Cell::setBField(double &B, double &theta, double &phi)
{
  //B field wrt to cell axes
  BField = B;
  double r=1;
  (BVectSpherical)[0]=r;
  (BVectSpherical)[1]=theta;
  (BVectSpherical)[2]=phi;
  //std::cout<<r<<"\t"<<theta<<"\t"<<phi<<"\n";
  
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Cell::setNeighbour(unsigned &side, Cell *c)
{ 
  if (side < 6)
    {
      // std::cout<<c->cellPos()<<"\t"<<this->cellPos()<<std::endl;
      
      neighbours[side] = c;
    }
  else
    {
      std::cerr<<"Cell Out of Range"<<std::endl;
    }
  
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Cell::outputNeighbours()
{
  std::cout<<"\n"<<this->cellPos()<<"\t"<<"cell"<<std::endl;
  Cell *temp1, *temp2, *temp3, *temp4, *temp5, *temp6;
  temp1 = neighbours.hexFront();
  if(temp1 != NULL)
    { std::cout<<temp1->cellPos()<<" neighbour"<<"\n";
    }
  temp2 = neighbours.hexBack();
  if(temp2 != NULL)
    { std::cout<<temp2->cellPos()<<" neighbour"<<"\n";
    }
  temp3 = neighbours.hexTop();
  if(temp3 != NULL)
    { std::cout<<temp3->cellPos()<<" neighbour"<<"\n";
    }
  temp4 = neighbours.hexBottom();
  if(temp4 != NULL)
    { std::cout<<temp4->cellPos()<<" neighbour"<<"\n";
    }
  temp5 = neighbours.hexRight();
  if(temp5 != NULL)
    { std::cout<<temp5->cellPos()<<" neighbour"<<"\n";
    }
  temp6 = neighbours.hexLeft();
  if(temp6 != NULL)
    { std::cout<<temp6->cellPos()<<" neighbour"<<"\n";
    }
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
Cell * Cell::retNeighbour(int &face)
{
  Cell *vicinus;
  
  if(face == 0)
    {
      vicinus = neighbours.hexFront();
    }
  else if(face == 1)
    {
      vicinus = neighbours.hexBack();
    }
  else if(face == 2)
    {
      vicinus = neighbours.hexTop();
    
    }
  else if(face == 3)
    {
     vicinus = neighbours.hexBottom();
   
    }
  else if(face == 4)
    {
      vicinus = neighbours.hexRight();
   
    }
  else if(face == 5)
    {
      vicinus = neighbours.hexLeft();
    }
  return vicinus;
  
}
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Cell::setPhotDirections()
{
  int g, alp, thet, face;
  double al, th;
  double r = 1;
  grid = pho->gridSize();
  agrid = pho->agridSize();
  
  for (g = 0; g < grid; ++g) 
    {
      for (alp = 0; alp < agrid; ++alp) 
	{
	  al = pho->phiElement(alp);
	  
	  for (thet = 0; thet < agrid; ++thet)
	    {
	      th = pho->thetaElement(thet);
	      face = whichFace(r, th, al);
	      pho->writeDirection(g,alp,thet,face);     
	    }
	}
    }
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
int Cell::whichFace(double &r, double &phi, double &theta)
{
  //pitch angle phi wrt Bz
  double xB, yB, zB, x, y, z;
  nums.spheriToCarte(r,phi,theta,x,y,z);
  int face = faceIntersect(x, y, z);
  return face;
  
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
int Cell::faceIntersect(double &x, double &y, double &z)
{
  //plane satisfies n.(r-r0) = 0, where n is normal, r and r0 two points on a plane
  //the vector from r to r0 will be normal to n, if r, r0 on the plane

  /* Forward face plane n = [1,0,0], r0 = [10,0,0]
     Back n = [-1,0,0], r0 = [-10,0,0]
     Left n = [0,0,-1], r0 = [0,0,-10]
     Right n = [0,0,1], r0 = [0,0,10]
     Up n = [0,1,0], r0 = [0,10,0]
     Down n = [0,-1,0], r0 = [0,-10,0] */

  /* Line through point r1 and r2. r1 will be from the simulation and r2 = r1*10 
   to put the face plane between r1 and r2.  Line is defined by 
   r = r1 + u(r2 -r1), where u is real. 
   Thus the plane and the line will intersect when:
   n.(r1 + u(r2 - r1)) = n.r0  

   Therefore u = (n.(r0-r1))/(n.(r2-r1))
   if the denominator is > 0 then the lines interesect
   if 0 < u < 1, then the intersection is between r1 and r2 */
 
  /*Note: some directions may intersect more than one orthogonal face planes, therefore
   need to check if the intersection is between two points.
  Or: With u available, the point of intersection r = r1 + u(r2 -r1), which can then be 
  checked to see which face plane it lies on using n.(r-r0) = 0 */

  //the following is for the exceptional cases when the photon/electrons are pointed to 
  //the edge of a cube. To resolve this issue, using a random number to fractionally
  //change the direction, therefore having a definite cube face for direction.
  int idum= (-1*time(NULL));
  if(x == 1.0  || x == 0.0 || x == -1.0)
    {
      x-=((nums.random(idum))/1000.0);
    }
  if(y == 1.0  || y == 0.0 || y == -1.0)
    {
      y-=((nums.random(idum))/1000.0);
    }
  if(z == 1.0  || z == 0.0 || z == -1.0)
    {
      z-=((nums.random(idum))/1000.0);
    }
    
  double uF=0.0, uB=0.0, uL=0.0, uR=0.0, uU=0.0, uD=0.0, fd=0.0, bd=0.0, 
    ld=0.0, rd=0.0, ud=0.0, dd=0.0;
  double r2x=10.*x, r2y=10.*y, r2z=10.*z;
  double rArr [] = {0.0,0.0,0.0};
  
  //set the cell boundaries at 10 and -10 = arbitrary scaling
  //forward plane, x =1 y=z=0 for the normal vector
  fd = (1.*(r2x-x));//n.r2-r1 => [1,0,0].([r2x,r2y,r2z]-[r1x,r1y,r1z])
                     //=> 1(r2x-r1x)+0(r2y-r1y)+0(r2z-r1z) => 1.*(r2x-r1x)
  if (fd > 0.0)
    {
      //r0 is point on a plane
      //on the forward plance (10,0,0) is on the plane
      //r0-r1 => [(10-x),(0-y),(0-z)]
      //n.(r0-r1)=>[1,0,0].[10-x,-y,-z]=>1.*(10-x)+0*(-y)+0*(-z)=>1.*(10-x)
      uF = (1.*(10.-x))/fd;//distance from r1 to intersectin point
      rArr[0] = x+(uF*(r2x-x));//intersection point.
      rArr[1] = y+(uF*(r2y-y));
      rArr[2] = z+(uF*(r2z-z));
      if(rArr[1] < 10.0 && rArr[1] > -10.0 && rArr[2] < 10.0 && rArr[2] > -10.0)
	{
	  return 0;
	}
    }

  //back, x = -1, y=z=0
  bd = (-1.*(r2x-x));
  if (bd > 0.0)
    {
      uB = (-1.*(10.-x))/bd;
      rArr[0] = x+(uB*(r2x-x));
      rArr[1] = y+(uB*(r2y-y));
      rArr[2] = z+(uB*(r2z-z));
      if(rArr[1] < 10.0 && rArr[1] > -10.0 && rArr[2] < 10.0 && rArr[2] > -10.0)
	{
	  return 1;
	}
      
    }

  //up, y = 1, x=z=0
  ud = (1.*(r2y-y));
  if (ud > 0.0)
    {
      uU = (1.*(10.-y))/ud;
      rArr[0] = x+(uU*(r2x-x));
      rArr[1] = y+(uU*(r2y-y));
      rArr[2] = z+(uU*(r2z-z));
      if(rArr[0] < 10.0 && rArr[0] > -10.0 && rArr[2] < 10.0 && rArr[2] > -10.0)
	{
	  return 2;
	}  
    }
  
  //down, y = -1, x=z=0
  dd = (-1.*(r2y-y));
  if (dd > 0.0)
    {
      uD = (-1.*(10.-y))/dd;
      rArr[0] = x+(uD*(r2x-x));
      rArr[1] = y+(uD*(r2y-y));
      rArr[2] = z+(uD*(r2z-z));
      if(rArr[0] < 10.0 && rArr[0] > -10.0 && rArr[2] < 10.0 && rArr[2] > -10.0)
	{
	  return 3;
	}
    }
 
  //right, z=1, y=x=0
  rd = (1.*(r2z-z));
  if (rd > 0.0)
    {
      uR = (1.*(10.-z))/rd;
      rArr[0] = x+(uR*(r2x-x));
      rArr[1] = y+(uR*(r2y-y));
      rArr[2] = z+(uR*(r2z-z));
      if(rArr[0] < 10.0 && rArr[0] > -10.0 && rArr[1] < 10.0 && rArr[1] > -10.0)
	{
	  return 4;
	}
    }
 
  //left, z =-1, y=x= 0 
  ld = (-1.*(r2z-z));
  if (ld > 0.0)
    {
      uL = (-1.*(10.-z))/ld;
      rArr[0] = x+(uL*(r2x-x));
      rArr[1] = y+(uL*(r2y-y));
      rArr[2] = z+(uL*(r2z-z));
      if(rArr[0] < 10.0 && rArr[0] > -10.0 && rArr[1] < 10.0 && rArr[1] > -10.0)
	{
	  return 5;
	}
    }
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Cell::physicalProcesses()
{
  //synch = new Synchrotron();
  std::string app1 = "bef";
  std::string app = "aft";
  
  const double dx=length;
  
  synch->synchPhiTheta(ele, pho, BField, BVectSpherical, length, synchArray);
  if(ext_compton_bool)
    {
      compt->comptonEmission(ele, pho, length);
    }
  else
    {
      //Compton passes intrinsic to transiting photon distributions
      //if compton is off then have to do this manually
      passIntrPhToTrans(pho);
    }
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Cell::passIntrPhToTrans(Photons *pho)
{
  int itNu, itPhi, itTheta, face, intrinsic = 6;
  int gridS = pho->gridSize();
  int agridS = pho->agridSize();

  double intrINu;
  
  for (itNu = 0; itNu <= gridS-1; ++itNu)
    {
      for (itPhi = 0; itPhi <= agridS-1; ++itPhi)
	{
	  for (itTheta = 0; itTheta <= agridS-1; ++itTheta)
	    {
	      intrINu = pho->get3DdistElement(itNu, itPhi, itTheta, intrinsic);
	      //distribution to which the photons will be passed
	      face = pho->getDirection(itNu, itPhi, itTheta);
	      pho->writeArr3Trans(itNu, itPhi, itTheta, intrINu, face);
	     
	    }
	}
    }
  
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
Photons * Cell::retPhot()
{
  return pho;
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
Electrons * Cell::retElec()
{
  return ele;
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
Arr3D<int, double> * Cell::neighbourDist(int &face)
{
  Photons *phots;
  Cell *vicinus;
  
  int i = 6, output=7;
  if(face == 0)
    {
      vicinus = neighbours.hexFront();
      if(vicinus != NULL)
	{
	  phots = vicinus->retPhot();
	  return phots->transiting(i);
	}
      else
	{
	  phots = this->retPhot();
	  return phots->transiting(output);
	}
    }
  if(face == 1)
    {
      
      vicinus = neighbours.hexBack();
      if(vicinus != NULL)
	{
	  phots = vicinus->retPhot();
	  return phots->transiting(i);
	}
      else
	{
	  phots = this->retPhot();
	  return phots->transiting(output);
	}
    }
  if(face == 2)
    {
      vicinus = neighbours.hexTop();
      if(vicinus != NULL)
	{
	  phots = vicinus->retPhot();
	  return phots->transiting(i);
	}
      else
	{
	  phots = this->retPhot();
	  return phots->transiting(output);
	}
    }
  if(face == 3)
    {
     vicinus = neighbours.hexBottom();
     if(vicinus != NULL)
       {
	 phots = vicinus->retPhot();
	 return phots->transiting(i);
       }
     else
       {
	 phots = this->retPhot();
	 return phots->transiting(output);
       } 
    }
  if(face == 4)
    {
      vicinus = neighbours.hexRight();
      if(vicinus != NULL)
	{
	  phots = vicinus->retPhot();
	  return phots->transiting(i);
	}
      else
	{
	  phots = this->retPhot();
	  return phots->transiting(output);
	} 
    }
  if(face == 5)
    {
      vicinus = neighbours.hexLeft();
      if(vicinus != NULL)
	{
	  phots = vicinus->retPhot();
	  return phots->transiting(i);
	}
      else
	{
	  phots = this->retPhot();
	  return phots->transiting(output);
	} 
    }
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

void Cell::passPhotonsToIntr()
{
  int intrinsic = 6;
  Photons *ph = (this)->retPhot();
  Cell *neighbour;
 
  for(int i = 0; i < 6; ++i)
    {
      neighbour = this->retNeighbour(i);
      if(neighbour != NULL)
	{
	  this->passArrs(ph, i);
	}
      else
	{
	  ph->addArray(this->neighbourDist(i), ph->transiting(i));
	}
            
    }
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Cell::passArrs(Photons *ph, int &f)
{
  
  Cell *neighbour = this->retNeighbour(f);
  Arr3D<int, double> *fromDist = ph->transiting(f);

  
  int intrinsic = 6;
  Photons *neighPhot = neighbour->retPhot();
  std::vector<double> *neighPhi = neighPhot->retPhi();
  std::vector<double> *neighTheta = neighPhot->retTheta();
  Arr3D<int, double> *toDist = neighPhot->transiting(intrinsic);
  int gridS = ph->gridSize();
  int agridS = pho->agridSize();
  int i, j, k, phiIndex, thetaIndex;
  double phi, theta, nu, val, r=1.;
  double x, y, z, nPhi, nTheta, nR;
  
   
  for (i = 0; i < gridS; ++i)
    {
      for (j = 0; j < agridS; ++j)
	{
	  phi = ph->phiElement(j);
	  
	  for (k = 0; k < agridS; ++k)
	    {
	      theta = ph->thetaElement(k);
	      val = fromDist->getElement(i, j, k);
	      phiIndex = getind.result(neighPhi, phi);
	      thetaIndex = getind.result(neighTheta, theta);
	      toDist->addElement(i, phiIndex, thetaIndex, val);
	    }
	}
    }
}


//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

void Cell::createSynchArray()
{
  int gridS = pho->gridSize();
  int agridS = pho->agridSize();
  synchArray = new Arr4D<int, double>(gridS, agridS, agridS, gridS, 0.0);
  int n, al, th, g;
  double alpha, phi, theta, nu, gamma, nuCrit,  x, val, sineAlpha;
  std::vector<double>::const_iterator nuIt, gIt;
  double con2 = (3*physcon.e*BField)/(4.*physcon.pi*physcon.m_e);//in SI units

  for (n = 0, nuIt = pho->xBegin(); n < gridS;
       ++n, ++nuIt) 
    {
      nu = *nuIt;
      
      for (al = 0; al < agridS; ++al) 
	{
	  phi = pho->phiElement(al);
	  for (th = 0; th < agridS; ++th)
	    {
	      theta = pho->thetaElement(th);
	      
	      for (g = 0, gIt = ele->xBegin(); g < gridS; ++g, ++gIt)
		{
		  gamma = *gIt;
		  alpha = nums.pitchAngle(phi, (BVectSpherical)[1], theta, (BVectSpherical)[2]);
		  sineAlpha = sin(alpha);
		  //sine become negative after pi. the following just repeats
		  //0 -> pi behaviour of the sine function
		  if(sineAlpha < 0.0) sineAlpha *= -1.;
		  nuCrit = con2 * sqr(gamma) * sineAlpha;
		  
		  if(nuCrit == 0.0)
		    {
		      x=0.0;
		    }
		  else
		    {
		      x = nu/nuCrit;		  
		    }
		  val = x * synch->synRad(x);
		  synchArray->putElement(n, al, th, g, val);
		}
	    }
	}
    }
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Cell::dopplerFactor()
{
  double lorentzF = ext_lorentz_factor;
  double beta = sqrt(1. - (1./sqr(lorentzF)));
  double viewTheta = ext_source_angle * 0.0174532925;//convert to radians
  double dopplerF = 1./(lorentzF * (1-beta*(cos(viewTheta))));
  return dopplerF;
  
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
int Cell::emptyFaceCount()
{
  int count = 0;
  Cell *temp1, *temp2, *temp3, *temp4, *temp5, *temp6;
  temp1 = neighbours.hexFront();
  if(temp1 == NULL)
    { 
      ++count;
    }
  temp2 = neighbours.hexBack();
  if(temp2 == NULL)
    { 
      ++count;
    }
  temp3 = neighbours.hexTop();
  if(temp3 == NULL)
    { 
      ++count;      
    }
  temp4 = neighbours.hexBottom();
  if(temp4 == NULL)
    { 
      ++count;
    }
  temp5 = neighbours.hexRight();
  if(temp5 == NULL)
    { 
      ++count;
    }
  temp6 = neighbours.hexLeft();
  if(temp6 == NULL)
    { 
      ++count;
    }
  return count;
  
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Cell::angleAberration(double thetaPrime)
{
  double lorentzF = ext_lorentz_factor;
  double beta = sqrt(1. - (1./sqr(lorentzF)));
  double cosThetaPrime = cos(thetaPrime);
  double cosTheta = (cosThetaPrime + beta)/(1.+(beta*cosThetaPrime));
  double val = acos(cosTheta);
  if(thetaPrime < 0.0)
    {
      val *= -1.0;
    }
  return val;
}
