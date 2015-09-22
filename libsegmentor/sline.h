// sline.hxx


#ifndef LIBSEGMENTOR_SLINE
#define LIBSEGMENTOR_SLINE

#include "dualnumber.h"
#include "ndm.h"

#include <iostream>
#include <math.h>

class sline{
private:
  dual sl[3];
public:
  sline ();
  sline(const double *pnt, const double *dir, int pflag=0);
  // Mostly this constructor is used. Just a point and a dir
  // or two points.
  // If the pnt is NULL then a line in the infinity.
  // If pflag != 0 then both are points on the line, if it is 0
  // then the second vector is the direction of the line.
  //  sline(double rho, double phi, double theta, double alpha);
  sline(const dual *sl);
  sline(dual x, dual y, dual z);
  ~sline();
  void reverse();
  // Negates direction
  const dual & component(int i) const;
  double direction(int i) const;
  void direction(double *d) const;
  double moment(int i) const;
  void moment(double *d) const;
  void nearp(double *p) const;
  double dist(const double *refp) const;
  double polar(double *plr) const;
  // It is not the same as parameterizing the cylinder, but it is its "dual"
  // The direction of the nearest point is expressed in terms of two
  // directions being orthogonal to the line.

  int finite() const;


friend dual      operator  % (const sline &a, const sline &b);
  // Dual angle of two straight lines = cos(f) - eps*g*sin(f), where
  // f is the real angle and g is the length of the transversal (!!!!)
  // The number:
  // cos(f) - eps*g*sin(f) = cos (f + eps*g)
  // by definition
  //
}; 

inline sline::sline()
{
}

inline int sline::finite() const
{
  return (sl[0].Dual() != HUGE);
}

std::ostream & operator << (std::ostream &os, const sline &s);

sline     operator  * (const sline &a, const sline &b);
// Transversal line

int rotax(const sline &s0,
	  const sline &s1,
	  const sline &s2,
	  const sline &s3,
	  sline *ressl,
	  double curvradmax = 1.0e+8,
	  double mindist    = 1.0e-3,
	  double planeps    = 1.0e-4);


int guess_axes(
 const sline *norms,      // Array of normal lines of a surface
 int          nnorms,     // Number of normals > 4
 sline      **axs,        // Output sline array address will be put here by
 // the  function. The order is determined by the number of cardinality in a
 // cell, the standard deviation of the 'nearest points to origin' and by
 // the standard deviation of the axis directions. The length of this array
 // is returned by the 'guess_axes' fuction. The array for output is
 // normally taken from the heap by 'new'. HOWEVER:
 // *!*
 // If below 'axsize' is not NULL and *axsize > 0, then it is assumed
 // that there is an array beginning at '*axs' which may contain 'axsize'
 // piece of 'sline'-s. In this case AT MOST '*axs' pieces are output and
 // copied to '*axs' and NO new is performed.
 // Otherwise if at the beginning it has been '*axsize <= 0' then the array
 // is taken from the heap by 'new' its address will written to 'axs', its
 // length will be returned by this function, AND this length (in terms of
 // 'sline'-s) will be written onto 'axszize' as well!!
 //
 int         *axsize=NULL, // The size info for 'axs' above.
 double       maxsdn = 10.0,   // Maximum st. deviation in "near" positions
 double     **sdns=NULL,  // Standard deviation of "near" positions (output)
 int          kslice=10,  // Number of cell slices on a halfsphere
                          // k*(k+1)/2 cells in total
 double       curvradmax = 1.0e+8, // Maximum values for t0, t1
                                   // only for "finite" axes
 double       mindist    = 1.0e-3, // Minimum meaningful distance
 double       planeps    = 1.0e-4, // Minimum value for angles or such
 int          mincno=2,   // Min. number of cardinality in a cell
 int        **cnos=NULL,  // Cardinality of cells (output)
 double       maxsdd =0.5,  // Maximum st. deviation in direction
                          // as a proportion to 3.545/k (2sqrt(M_PI)/k)
 double     **sdds=NULL); // St. dev of directions (output)
//
// Returns the number of axis which have met all the criteria
// in the order of the cardinality within the cell. The area where
// the solutions are put is taken from the heap by 'new'.
// the same applies for information arrays.
//
#endif

