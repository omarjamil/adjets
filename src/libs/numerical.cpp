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
#include <limits>

#include "numerical.hpp"

#define sqr(a) ((a)*(a))
//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

Numerical::Numerical()
{}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

Numerical::~Numerical()
{}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
//Chebeyshev Approximation
//Numerical recipes page 199
double Numerical::chebev(const double a, const double b, std::vector<double> &c, 
	    const int m, const double x) 
{
  double d=0.0, dd=0.0, sv, y,y2;
  int j;
  
  if ((x-a)*(x-b) > 0.0)
    {
      std::cerr<<"x not in range in routine chebev"<<std::endl;
    }
  
  y2=2.0*(y=(2.0*x-a-b)/(b-a));
  
  for (j=m-1;j>=1;--j) 
    { 
      sv=d; 
      d=y2*d-dd+c[j]; 
      dd=sv; 
    }
  return y*d-dd+0.5*c[0];
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
//NR page 250
void Numerical::beschb(double x, double *gam1, double *gam2, 
	    double *gampl, double *gammi)
{
  const int NUSE1=7, NUSE2=8;
  static const double c1_d[7] = {
    -1.142022680371168e0,6.5165112670737e-3, 3.087090173086e-4,
    -3.4706269649e-6,6.9437664e-9, 3.67795e-11,-1.356e-13};

  static const double c2_d[8] = { 
    1.843740587300905e0,-7.68528408447867e-2, 1.2719271366546e-3,
    -4.9717367042e-6,-3.31261198e-8, 2.423096e-10,-1.702e-13, 
    -1.49e-15};
  double xx;
  static std::vector<double> c1(c1_d, c1_d+7);
  static std::vector<double> c2(c2_d, c2_d+8);
  
  xx=8.0*x*x-1.0;
  *gam1=chebev(-1.0,1.0,c1,NUSE1,xx);
  *gam2=chebev(-1.0,1.0,c2,NUSE2,xx);
  *gampl= *gam2-x*(*gam1);
  *gammi= *gam2+x*(*gam1);
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....

//Modified Bessel function NR page 253
// void Numerical::bessik(const double x, const double xnu, 
// 		       double &ri, double &rk, double &rip, 
// 		       double &rkp)

void Numerical::bessik(const double x, const double xnu, 
		       double &rk)
{
  const int MAXIT=10000;
  const double EPS=std::numeric_limits<double>::epsilon();
  const double FPMIN=std::numeric_limits<double>::min()/EPS;
  const double XMIN=2.0, PI=3.141592653589793;
  
  double a, a1, b, c, d, del, del1, delh, dels, e, f, fact, fact2, 
    ff, gam1, gam2, gammi, gampl, h, p, pimu, q, q1, q2, qnew, ril, 
    ril1, rimu, rip1, ripl, ritemp, rk1, rkmu, rkmup, rktemp, s, sum, sum1,
    x2, xi, xi2, xmu, xmu2;
  
  int i, l, nl;
  
  if (x <= 0.0 || xnu < 0.0) 
    {
      std::cerr<<"bad arguments in bessik"<<std::endl;
    }
  
  nl=int(xnu+0.5);
  xmu=xnu-nl;
  xmu2=xmu*xmu;
  xi=1.0/x;
  xi2=2.0*xi;
  h=xnu*xi;
  
  if (h < FPMIN)
    {
      h=FPMIN;
    }
  
  b=xi2*xnu;
  d=0.0;
  c=h;
  d=0.0;
  c=h;
  
  for (i=1;i<=MAXIT;i++) 
    { 
      b += xi2; 
      d=1.0/(b+d); 
      c=b+1.0/c; 
      del=c*d; 
      h=del*h; 
      if (fabs(del-1.0) < EPS) break; 
    }
  
  if (i > MAXIT) //"x too large in bessik; try asymptotic expansion"
    {
      std::cerr<<"x too large in bessik; try asymptotic expansion"
	       <<std::endl;
    }
  
  ril=FPMIN;
  ripl=h*ril;
  ril1=ril;
  rip1=ripl;
  fact=xnu*xi;

  for (l=nl;l>=1;l--) 
    { 
      ritemp=fact*ril+ripl; 
      fact -= xi;
      ripl=fact*ritemp+ril; 
      ril=ritemp; 
    }
  
  f=ripl/ril;
  
  if (x < XMIN) 
    {
    x2=0.5*x;
    pimu=PI*xmu;
    fact = (fabs(pimu) < EPS ? 1.0 : pimu/sin(pimu));
    d = -log(x2);
    e=xmu*d;
    fact2 = (fabs(e) < EPS ? 1.0 : sinh(e)/e);
    beschb(xmu,&gam1,&gam2,&gampl,&gammi);
    ff=fact*(gam1*cosh(e)+gam2*fact2*d);
    sum=ff;
    e=exp(e);
    p=0.5*e/gampl;
    q=0.5/(e*gammi);
    c=1.0;
    d=x2*x2;
    sum1=p;
    
    for (i=1;i<=MAXIT;i++) 
      {
	ff=(i*ff+p+q)/(i*i-xmu2);
	c *= (d/i);
	p /= (i-xmu);
	q /= (i+xmu);
	del=c*ff;
	sum += del;
	del1=c*(p-i*ff);
	sum1 += del1;
	
	if (fabs(del) < fabs(sum)*EPS)
	  {
	    break;
	  }
      }
    
    if (i > MAXIT)
      {
	std::cerr<<"bessk series failed to converge"<<std::endl;
      }
    rkmu=sum;
    rk1=sum1*xi2;
    }
  else 
    {
      b=2.0*(1.0+x);
      d=1.0/b;
      h=delh=d;
      q1=0.0;
      q2=1.0;
      a1=0.25-xmu2;
      q=c=a1;
      a = -a1;
      s=1.0+q*delh;
    
      for (i=2;i<=MAXIT;i++) 
	{
	  a -= 2*(i-1);
	  c = -a*c/i;
	  qnew=(q1-b*q2)/a;
	  q1=q2;
	  q2=qnew;
	  q += c*qnew;
	  b += 2.0;
	  d=1.0/(b+a*d);
	  delh=(b*d-1.0)*delh;
	  h += delh;
	  dels=q*delh;
	  s += dels;
	  
	  if (fabs(dels/s) < EPS) 
	    {
	      break;
	    }
	}
      
      if (i > MAXIT) 
	{
	  std::cerr<<"bessik: failure to converge in cf2"<<std::endl;
	}
      
      h=a1*h;
      rkmu=sqrt(PI/(2.0*x))*exp(-x)/s;
      rk1=rkmu*(xmu+x+0.5-h)*xi;
    }
  
  // rkmup=xmu*xi*rkmu-rk1;
  // rimu=xi/(f*rkmu-rkmup);
  // ri=(rimu*ril1)/ril;
  // rip=(rimu*rip1)/ril;
  
  for (i=1;i<=nl;i++) 
    { 
      rktemp=(xmu+i)*xi2*rk1+rkmu; 
      rkmu=rk1; 
      rk1=rktemp; 
    }
  
  rk=rkmu;
  // rkp=xnu*xi*rkmu-rk1;
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Numerical::gammalnFunction(const double &xx)
{
  int j;
  double x,tmp,y,ser;
  static const double cof[14]={57.1562356658629235,-59.5979603554754912,
			     14.1360979747417471,-0.491913816097620199,.339946499848118887e-4,
			     .465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3,
			     -.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,
			     .844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5};
  if (xx <= 0) 
    {
      std::cerr<<"bad argument in Numerical::gammalnFunction"<<std::endl;;
    }
  
  y=x=xx;
  tmp = x+5.24218750000000000;
  tmp = (x+0.5)*log(tmp)-tmp;
  ser = 0.999999999999997092;
  for (j=0;j<14;j++) 
    {
      ser += cof[j]/++y;
    }
  return tmp+log(2.5066282746310005*ser/x);
 
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
//From Numerical Recipes ran2
double Numerical::random(int &idum)
{
  const int im1=2147483563, im2=2147483399;
  const int ia1=40014, ia2=40692, iq1=53668, iq2=52774, ir1=12211;
  const int ir2=3791, ntab=32, imm1=(im1-1);
  const int ndiv=(1+imm1/ntab);
  const double eps=1.2e-7, RNMX=(1.0-eps), am=(1./double(im1));
  int j, k;
  static int idum2=123456789;
  static int iy=0;
  static std::vector<int> iv(ntab);
  double temp;
  
  if (idum <= 0) {
    if (-(idum) < 1) idum=1;
    else idum = -(idum);
    idum2=(idum);
    for (j=ntab+7;j>=0;--j) {
      k=(idum)/iq1;
      idum=ia1*(idum-k*iq1)-k*ir1;
      if (idum < 0) idum += im1;
      if (j < ntab) iv[j] = idum;
    }
    iy=iv[0];
  }
  k=(idum)/iq1;
  idum=ia1*(idum-k*iq1)-k*ir1;
  if (idum < 0) idum += im1;
  k=idum2/iq2;
  idum2=ia2*(idum2-k*iq2)-k*ir2;
  if (idum2 < 0) idum2 += im2;
  j=iy/ndiv;
  iy=iv[j]-idum2;
  iv[j] = idum;
  if (iy < 1) iy += imm1;
  if ((temp=am*iy) > RNMX) return RNMX;
  else return temp;
  
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Numerical::spheriToCarte(double &r, double &theta, double &phi, 
			      double &x, double &y, double &z)
{
  // y = r*cos(theta)*sin(phi);
  // z = r*sin(theta)*sin(phi);
  // x = r*cos(phi);

  x = r*sin(theta)*sin(phi);
  y = r*cos(phi);
  z = r*cos(theta)*sin(phi);
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
void Numerical::carteToSpheri(double &x, double &y, double &z,
			      double &r, double &theta, double &phi)
{
  r = sqrt(sqr(x)+sqr(y)+sqr(z));
  int quad = whichQuadrant(x,y,z);
  
  if(quad >= 1 && quad <= 4)
    {
      if(z == 0.) 
	{
	  theta = 1.57079632679;
	}
      else 
	{
	  theta = atan(x/z);
	  if(theta < 0.0) theta += 3.14159265358;
	}
      phi = acos(y/r);
      
      if(phi < 0.0) phi *= -1.;
    }
  else if(quad >= 5 && quad <= 8)
    {
      if(z == 0.) 
	{
	  theta = 1.5 * 3.14159265358;
	}
      else 
	{
	  theta = atan(x/z);
	  if(theta > 0.0) theta += 3.14159265358;
	  if(theta < 0.0) theta += (2.*3.14159265358);
	       
	}
      phi = asin(y/r) + (1.5*3.14159265358);
      if(phi < 0.0) phi *= -1.;
    }
  else 
    {
      std::cerr<<"Error in carteToSpheri \n";
    }
   
}


//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
int Numerical::whichQuadrant(double &x, double &y, double &z)
{
  //see notebook for quadrant geomtery, but it is basically 
  //starting from positive x, staying there going clockwise.
  //then to negative x, but positive y and z, and once again
  //going clockwise.
  int quad=0;
  if(x >= 0.0 && y >= 0.0 && z >= 0.0)
    {
      quad = 1;
    }
  else if(x >= 0.0 && y >= 0.0 && z < 0.0)
    {
      quad = 2;
    }
  else if(x >= 0.0 && y < 0.0 && z >= 0.0)
    {
      quad = 3;
    }
  else if(x >= 0.0 && y < 0.0 && z < 0.0)
    {
      quad = 4;
    }
  else if(x < 0.0 && y >= 0.0 && z >= 0.0)
    {
      quad = 5;
    }
  else if(x < 0.0 && y >= 0.0 && z < 0.0)
    {
      quad = 6;
    }
  else if(x < 0.0 && y < 0.0 && z >= 0.0)
    {
      quad = 7;
    }
  else if(x < 0.0 && y < 0.0 && z <= 0.0)
    {
      quad = 8;
    }
    else 
      {
	std::cerr<<"Error in which Quadrant \t"<<x<<"\t"<<y<<"\t"<<z<<"\n";
      }
    
    return quad;
    
}

//.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....oooOO0OOooo.....
double Numerical::pitchAngle(double &phi, double &phiB, double &theta, 
		  double &thetaB)
{
  double angleC = thetaB-theta;
  if(angleC < 0.0)
    {
      angleC = (2.*3.14159265358) + angleC;
    }
  double cosAlpha = (cos(phi)*cos(phiB)) + (sin(phi)*sin(phiB)*cos(angleC));
  return acos(cosAlpha);
  
}
