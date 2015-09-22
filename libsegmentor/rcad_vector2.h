/* rcad_vector.h */

#ifndef LIBSEGMENTOR_VECTOR_H
#define LIBSEGMENTOR_VECTOR_H 1

#ifdef __cplusplus
extern "C" {
#endif

#include <float.h>
#include <math.h>


/*---------------------------------------------------------------------------*/

inline double
s_mult(const double *a, const double *b)
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

/*---------------------------------------------------------------------------*/

inline void
assign( double *a, const double *b)
{
  a[0] = b[0];
  a[1] = b[1];
  a[2] = b[2];
}

/*---------------------------------------------------------------------------*/

inline void
v_mult( double *r, const double *a, const double *b)
{
  double temp[3];

  temp[0] = a[1]*b[2] - a[2]*b[1];
  temp[1] = a[2]*b[0] - a[0]*b[2];
  temp[2] = a[0]*b[1] - a[1]*b[0];

  assign(r,temp); 

}

/*---------------------------------------------------------------------------*/

inline void
mult( double *r, double s, const double *v)
{
  r[0] = s*v[0];
  r[1] = s*v[1];
  r[2] = s*v[2];
}

/*---------------------------------------------------------------------------*/

inline void
sub( double *r, const double *a, const double *b)
{
  r[0] = a[0]-b[0];
  r[1] = a[1]-b[1];
  r[2] = a[2]-b[2];
}

/*---------------------------------------------------------------------------*/
inline void
add( double *r, const double *a, const double *b)
{
  r[0] = a[0]+b[0];
  r[1] = a[1]+b[1];
  r[2] = a[2]+b[2];
}

/*--------------------------------------------------------------------------*/ 
inline double
vabs( double *r)
{
  return sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
}

/*--------------------------------------------------------------------------*/ 
inline void
unit( double *r)
{
  double  d;

  d = vabs( r);

  r[0] /= d;
  r[1] /= d;
  r[2] /= d;
}

/*--------------------------------------------------------------------------*/ 
inline static	void
makeCoordinateSystem( double *z, double *x, double *y)
{

   double	ax;
   double      ay;
   double      az;
   double	edir[3];

   edir[0] = edir[1] = edir[2] = 0.0;
   ax = fabs(z[0]);
   ay = fabs(z[1]);
   az = fabs(z[2]);

   if (ax > ay)
   {
      if (ax > az)
      {
	 if (ay > az)
	    edir[2] = 1.0;
	 else
	    edir[1] = 1.0;
      }
      else
	 edir[1] = 1.0;
   }
   else {
      if (ax < az)
	 edir[0] = 1.0;
      else
	 edir[2] = 1.0;
   }

   unit( z);
   v_mult( x, z, edir);
   unit( x);
   v_mult( y, z, x);
   unit( y); 
}
/*--------------------------------------------------------------------------*/ 
#ifdef __cplusplus
}
#endif

#endif





