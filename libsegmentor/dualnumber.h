// dualnumber.hxx


#ifndef LIBSEGMENTOR_DUALNUMBER
#define LIBSEGMENTOR_DUALNUMBER

#include <iostream>


class dual{
private:
  double      Re;      // Real part
  double      Du;      // Pure dual part
public:
  dual(double re = 0.0, double du=0.0);
  dual(const dual&);
  // ~dual();

  double Real() const;
  double Dual() const;

  dual& operator= (double);
  dual& operator+= (double);
  dual& operator-= (double);
  dual& operator*= (double);
  dual& operator/= (double);

  dual& operator= (const dual&);
  dual& operator+= (const dual&);
  dual& operator-= (const dual&);
  dual& operator*= (const dual&);
  dual& operator/= (const dual&);
}; 

dual operator+(const dual&, const dual&);
dual operator+(const dual&, const double&);
dual operator+(const double&, const dual&);


dual operator-(const dual&, const dual&);
dual operator-(const dual&, const double&);
dual operator-(const double&, const dual&);


dual operator*(const dual&, const dual&);
dual operator*(const dual&, const double&);
dual operator*(const double&, const dual&);


dual operator/(const dual&, const dual&);
dual operator/(const dual&, const double&);
dual operator/(const double&, const dual&);

dual operator+(const dual&);
dual operator-(const dual&);


dual sqrt (const dual& d);
dual cos (const dual&);
dual sin (const dual&);

dual acos (const dual&);
dual asin (const dual&);


const dual duone(1.0);
const dual duepsilon(0.0,1.0);


inline double dual::Real() const
{return Re;}

inline double dual::Dual() const
{return Du;}


inline dual::dual(double re, double du) :
Re(re),
Du(du)
{}

inline dual::dual(const dual& dd) :
Re(dd.Real()),
Du(dd.Dual())
{}

inline double Real(const dual& a)
{
  return(a.Real());
}

inline double Dual(const dual& a)
{
  return(a.Dual());
}

std::ostream & operator << (std::ostream &os, const dual &d);


dual s_mult(const dual *a, const dual *b);
void assign( dual *a, const dual *b);
void v_mult( dual *r, const dual *a, const dual *b);
void mult( dual *r, dual s, const dual *v);
void sub( dual *r, const dual *a, const dual *b);
void add( dual *r, const dual *a, const dual *b);
dual vabs( dual *r);
void unit( dual *r);

/*
double s_mult(const double *a, const double *b);
void assign( double *a, const double *b);
void v_mult( double *r, const double *a, const double *b);
void mult( double *r, double s, const double *v);
void sub( double *r, const double *a, const double *b);
void add( double *r, const double *a, const double *b);
double vabs( double *r);
void unit( double *r);
*/

#endif


