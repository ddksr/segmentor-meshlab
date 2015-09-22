//  dualnumber.cxx

#include "dualnumber.h"
#include "rcad_vector2.h"

#include <math.h>

// dual:: ~dual();


dual& dual::operator = (double re)
{
  this->Re = re;
  this->Du = 0.0;
  return *this;
}

dual& dual::operator+= (double re)
{
  this->Re += re;
  return *this;
}

dual& dual::operator-= (double re)
{
  this->Re -= re;
  return *this;
}


dual& dual::operator*= (double re)
{
  this->Re *= re;
  this->Du *= re;
  return *this;
}

dual& dual::operator/= (double re)
{
  this->Re /= re;
  this->Du /= re;
  return *this;
}

dual& dual::operator= (const dual& d)
{
  this->Re = d.Real();
  this->Du = d.Dual();
  return *this;
}


dual& dual::operator+= (const dual& d)
{
  this->Re += d.Real();
  this->Du += d.Dual();
  return *this;
}

dual& dual::operator-= (const dual& d)
{
  this->Re -= d.Real();
  this->Du -= d.Dual();
  return *this;
}


dual& dual::operator*= (const dual& d)
{
  this->Du = this->Re*d.Dual() + this->Du*d.Real();
  this->Re *= d.Real();
  return *this;
}

dual& dual::operator/= (const dual& d)
{
  this->Du = (this->Du*d.Real() - this->Re*d.Dual())
             / (d.Real()*d.Real());
  this->Re /= d.Real();
  return *this;
}

dual operator+(const dual& a, const dual& b)
{
  return dual(a.Real()+b.Real(),a.Dual()+b.Dual());
}

dual operator+(const dual& a, const double& b)
{
  return dual(a.Real()+b, a.Dual());
}

dual operator+(const double& b, const dual& a)
{
  return dual(a.Real()+b, a.Dual());
}

dual operator-(const dual&a, const dual&b)
{
  return dual(a.Real()-b.Real(),a.Dual()-b.Dual());
}

dual operator-(const dual& a, const double& b)
{
  return dual(a.Real()-b, a.Dual());
}

dual operator-(const double& b, const dual& a)
{
  return dual(b-a.Real(), -a.Dual());
}


dual operator*(const dual& a, const dual& b)
{
  return dual(a.Real()*b.Real(),
	      a.Real()*b.Dual() + a.Dual()*b.Real());
}

dual operator*(const dual& a, const double& b)
{
  return dual(a.Real()*b, a.Dual()*b);
}

dual operator*(const double& b, const dual& a)
{
  return dual(a.Real()*b, a.Dual()*b);
}


dual operator/(const dual& a, const dual& b)
{
  return dual(a.Real() / b.Real(),
	      (a.Dual()*b.Real() - a.Real()*b.Dual()) / (b.Real()*b.Real()));
}

dual operator/(const dual& a, const double& b)
{
  return dual(a.Real() / b,
	  a.Dual()/ b);
}

dual operator/(const double& a, const dual& b)
{
  return dual(a / b.Real(),
	      - a*b.Dual() / (b.Real()*b.Real()));
}

dual operator+(const dual& a)
{
  return dual(a.Real(),a.Dual());
}

dual operator-(const dual& a)
{
  return dual(-a.Real(),-a.Dual());
}

dual sqrt (const dual& d)
{
  double r = sqrt(d.Real());
  return dual( r,
	       r == 0.0 ? 1.0 : 0.5*d.Dual()/r);
}


dual cos (const dual& Phi)
{
  return dual( cos(Phi.Real()),
	       -Phi.Dual()*sin(Phi.Real()));
}

dual sin (const dual& Phi)
{
  return dual( sin(Phi.Real()),
	       Phi.Dual()*cos(Phi.Real()));
}

dual acos (const dual& d)
{
  double r=acos(d.Real());
  return dual( r,  d.Dual() == 0.0 ? 0.0 : -d.Dual()/sin(r));
}

dual asin (const dual& d)
{
  double r=asin(d.Real());
  return dual( r, d.Dual() == 0.0 ? 0.0 : d.Dual()/cos(r));
}


std::ostream & operator << (std::ostream &os, const dual &d)
{

  std::cout << '(' << d.Real() << ',' << d.Dual() << ')';

  return os;
}


// Vectorial operations


dual
s_mult(const dual *a, const dual *b)
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

/*---------------------------------------------------------------------------*/

void
assign( dual *a, const dual *b)
{
  a[0] = b[0];
  a[1] = b[1];
  a[2] = b[2];
}

/*---------------------------------------------------------------------------*/

void
v_mult( dual *r, const dual *a, const dual *b)
{
  dual temp[3];

  temp[0] = a[1]*b[2] - a[2]*b[1];
  temp[1] = a[2]*b[0] - a[0]*b[2];
  temp[2] = a[0]*b[1] - a[1]*b[0];

  assign(r,temp); 

}

/*---------------------------------------------------------------------------*/

void
mult( dual *r, dual s, const dual *v)
{
  r[0] = s*v[0];
  r[1] = s*v[1];
  r[2] = s*v[2];
}

/*---------------------------------------------------------------------------*/

void
mult( dual *r, double s, const dual *v)
{
  r[0] = s*v[0];
  r[1] = s*v[1];
  r[2] = s*v[2];
}

/*---------------------------------------------------------------------------*/

void
sub( dual *r, const dual *a, const dual *b)
{
  r[0] = a[0]-b[0];
  r[1] = a[1]-b[1];
  r[2] = a[2]-b[2];
}


/*---------------------------------------------------------------------------*/
void
add( dual *r, const dual *a, const dual *b)
{
  r[0] = a[0]+b[0];
  r[1] = a[1]+b[1];
  r[2] = a[2]+b[2];
}


/*--------------------------------------------------------------------------*/ 
dual
vabs( dual *r)
{
  return sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
}

/*--------------------------------------------------------------------------*/ 
void
unit( dual *r)
{
  dual  d;

  d = vabs( r);

  r[0] /= d;
  r[1] /= d;
  r[2] /= d;
}

// --------------------------------------------------------------------------

// Added by Bojan Kverh 25.2.1998
/*
double s_mult(const double *a, const double *b)
{ return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  }

void assign( double *a, const double *b) 
{
  a[0] = b[0];
  a[1] = b[1];
  a[2] = b[2];
}

void v_mult( double *r, const double *a, const double *b)
{

  r[0] = a[1]*b[2] - a[2]*b[1];
  r[1] = a[2]*b[0] - a[0]*b[2];
  r[2] = a[0]*b[1] - a[1]*b[0];

}

void mult( double *r, double s, const double *v)
{
  r[0] = s*v[0];
  r[1] = s*v[1];
  r[2] = s*v[2];
}

void sub( double *r, const double *a, const double *b)
{
  r[0] = a[0]-b[0];
  r[1] = a[1]-b[1];
  r[2] = a[2]-b[2];
}

void add( double *r, const double *a, const double *b)
{
  r[0] = a[0]+b[0];
  r[1] = a[1]+b[1];
  r[2] = a[2]+b[2];
}

double vabs( double *r)
{
  return sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
}

void unit( double *r)
{
  double  d;

  d = vabs( r);

  r[0] /= d;
  r[1] /= d;
  r[2] /= d;
}
*/
