// linegeo.cxx
// 9-May-97
// GL Cardiff
//

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include <iostream>
#include "rcad_vector2.h"
#include "dualnumber.h"
#include "linegeo.h"

#ifndef SMALLREAL
  #define SMALLREAL DBL_EPSILON
#endif

#define SMALL(a,b) (fabs(a) < b)

#define M_PI 3.14159265358979323846

#define DEBUG_TRACE

#ifdef DEBUG_TRACE
static int    debug_axes=0;
#endif

static const  double  degpi = 180.0/M_PI;

sline::sline(const double *pnt, const double *dir, int pflag)
{
  double mom[3];
  double p1p2[3];


  assign(p1p2,dir);
  if (pnt != NULL) {
    if (pflag !=0) {  // Both are points
      sub(p1p2, dir, pnt);
    }
    unit(p1p2);
    v_mult(mom, pnt, p1p2);

    sl[0] = dual(p1p2[0], mom[0]);
    sl[1] = dual(p1p2[1], mom[1]);
    sl[2] = dual(p1p2[2], mom[2]);
  } else {
    unit(p1p2);

    sl[0] = dual(p1p2[0], HUGE);
    sl[1] = dual(p1p2[1], HUGE);
    sl[2] = dual(p1p2[2], HUGE);
  }
}

//sline::sline(double rho, double phi, double theta, double alpha)
//{
//}


sline::sline(const dual *sli)
{
  int i;

  for (i=0; i < 3; ++i)
    sl[i] =sli[i];

  unit(sl);
}

sline::sline(dual x, dual y, dual z)
{
  sl[0] = x;
  sl[1] = y;
  sl[2] = z;
  unit(sl);
}

sline::~sline()
{}



const dual & sline::component(int i) const
{

  if ((i < 0) || (i >= 3)) {
    printf("Wrong indexing in sline::component:%d\n", i);
    i = 0;
  }
  return sl[i];
}

double sline::direction(int i) const
{
  if ((i < 0) || (i >= 3)) {
    printf("Wrong indexing in sline::direction:%d\n", i);
    i = 0;
  }
  return sl[i].Real();
}


void sline::direction(double *d) const
{
  d[0] = sl[0].Real();
  d[1] = sl[1].Real();
  d[2] = sl[2].Real();
}


double sline::moment(int i) const
{
  if ((i < 0) || (i >= 3)) {
    printf("Wrong indexing in sline::moment:%d\n", i);
    i = 0;
  }
  return sl[i].Dual();
}

void sline::moment(double *d) const
{
  d[0] = sl[0].Dual();
  d[1] = sl[1].Dual();
  d[2] = sl[2].Dual();
}


double sline::dist(const double *refp) const
{
  double refpxdir[3];
  double mom[3];

  direction(refpxdir);
  moment(mom);

  v_mult(refpxdir, refp, refpxdir);
  sub(refpxdir,refpxdir,mom);
  return vabs(refpxdir);
}
  


void sline::nearp(double *p) const
{
  double a[3];

  moment(p);

  if (finite()) {
    direction(a);
    v_mult(p,a,p);
  }
}


double sline::polar(double *plr) const
{
  double np[3];
  double dir[3];
  //  double rho,phi,theta,alpha;
  double bnf[3];
  double nt[3];
  double t,x,y;
  int i;
  // It is not the same as parameterizing the cylinder, but it is its "dual"
  // The direction of the nearest point is expressed in terms of two
  // directions being orthogonal to the line.

  direction(dir);
  plr[2] = acos(dir[2]);
  plr[1] = atan2(dir[1],dir[0]);

  if (finite()) {
    nearp(np);
    plr[0] = vabs(np);
    if (plr[0] > DBL_EPSILON) {
      mult(np, 1.0/plr[0], np);
    } else {
      t = fabs(dir[0]); i=0;
      if (fabs(dir[1]) < t) {
	t = fabs(dir[1]);
	i = 1;
      }
      if (fabs(dir[2]) < t) {
	i = 2;
      }
      np[0] = np[1] = np[2] = 0.0;
      np[i] = 1.0;
      v_mult(np, np, dir);
      unit(np);
      plr[0] = 0.0;
    }

    //  plr[2] = acos(np[2]);
    //  plr[1] = atan2(np[1],np[0]);

    bnf[0] = -sin(plr[1]);
    bnf[1] = cos(plr[1]);
    bnf[2] = 0.0;

    nt[0]  = cos(plr[2]);
    nt[1]  = -bnf[0]*nt[0];
    nt[0] *= bnf[1];
    nt[2]  = -sin(plr[2]);

  //  y = s_mult(bnf,dir);
  //  x = s_mult(nt, dir);

    y = s_mult(bnf,np);
    x = s_mult(nt, np);

    plr[3]=atan2(y,x);
  
  } else {
    plr[0] = HUGE;
    plr[3] = 0.0;
  }

  return plr[0];

}


dual      operator  % (const sline &a, const sline &b)
{
  // Dual angle of two straight lines = cos(f) - eps*g*sin(f), where
  // f is the real angle and g is the length of the transversal (!!!!)
  // The number:
  // cos(f) - eps*g*sin(f) = cos (f + eps*g)
  // by definition
  //
  return( s_mult(a.sl,b.sl));
}



sline     operator  * (const sline &a, const sline &b)
{
  // Transversal
  double adir[3];
  double bdir[3];
  double axbdir[3];
  double tdir[3];
  double pa[3];
  double pbpa[3];
  double *rpa;
  double ta;
  double nanb;

  a.direction(adir);
  b.direction(bdir);
  v_mult( axbdir, adir, bdir);
  
  if ((a.finite()) && (b.finite())) {
    a.nearp(pa);
    b.nearp(pbpa);
    sub(pbpa, pbpa, pa);

    nanb = s_mult(adir,bdir);
    mult (tdir, nanb, bdir);
    sub( tdir, adir, tdir);
  
    ta = s_mult(pbpa, tdir) / (1.0 - nanb*nanb);
    // Denominator should not be zero if adir is not parallel to bdir

    mult(adir, ta, adir);
    add (pa, pa, adir);
    rpa = pa;
  } else {
    rpa = NULL;
  }
  return sline(pa, axbdir);
}




std::ostream & operator << (std::ostream &os, const sline &s)
{
  double  plr[4];
  double  dir[3];
  double  near[3];

  s.polar(plr);
  s.direction(dir);
  s.nearp(near);

  std::cout << "///D: (" << dir[0] << ", " << dir[1] << ", " << dir[2] <<
    ") /NP:";
  if (s.finite()) {
    std::cout << near[0] << ", " << near[1] << ", " <<
      near[2] << ')' << std::endl;
    std::cout << "///rho:" << plr[0] << ", Phi:" << degpi*plr[1] << ", Theta:"
	 << degpi*plr[2] << ", Alpha: " << degpi*plr[3];
  } else {
    std::cout << "INFINITE" << std::endl;
    std::cout << "///rho:" << "INF" << ", Phi:" << degpi*plr[1] << ", Theta:"
	 << degpi*plr[2] << ", Alpha: " << degpi*plr[3];
  }

  return os;

}



static double triple(double *a, double *b, double *c)
{
  double axb[3];


  v_mult(axb, a, b);
  return (s_mult(axb,c));
}


static int        v_deg2eqn(double a,
		     double b,
		     double c,
		     double t1,
		     double t2,
		     double roottol,
		     //Root will be toleranced with this value.
		     // For near double roots
		     double *rp) // Results 
/* RETURN: Number of roots in the [t1,t2] interval. (-1) if there are
 * infinite number of roots.
 */
{
/* Second degree equation solver. Note that it uses the machine precision
 * 'SMALLREAL'. It solves the equation ax^2+bx+c=0 in the interval
 * [t1-roottol,t2+roottol].
 * COMPLEX roots are recognized if
 *      b^2/4a^2 - c/a < -max( roottol^2, 100*SMALLREAL)
 */

   double          discr,
                   x,
                   x2;
   double          eps_adiscr = 100.0 * SMALLREAL;	/* Magic number */
   double          eps_discr;
   int             rnum,
                   oldrnum;
   int             twook = 1;

   roottol = fabs(roottol);

   if ((eps_discr = roottol * roottol) < eps_adiscr)
      eps_adiscr = eps_discr;


   if (SMALL(a, SMALLREAL))
   {

      if (SMALL(b, SMALLREAL))
      {

	 if (SMALL(c, SMALLREAL))
	    return (-1);
	 else
	    return (0);

      }
      
      rp[0] = -c / b;
      rnum = 1;

   }
   else
   {

      x = -0.5 * (b / a);
      x2 = x * x;
      discr = x2 - c / a;		/* Hopefully does not loose accuracy */
      x2 *= 10.0 * SMALLREAL;

      if (eps_discr < x2)
	 eps_discr = x2;

      if ((discr < eps_adiscr) && (discr > -eps_discr))
      {
	 rp[0] = x;
	 rnum = 1;
      }
      else if (discr > 0.0)
      {

	 discr = sqrt(discr);
	 rp[0] = (fabs(x) < SMALLREAL ? x - discr :
		                        x + (x > 0.0 ? discr : -discr));
	 rp[1] = (fabs(rp[0]) > 1.0 ? c / (rp[0] * a) : 2.0 * x - rp[0]);
	 rnum = 2;

      }
      else
	 return (0);
   }

   oldrnum = rnum;

   if (rp[0] < t1 - roottol || rp[0] > t2 + roottol)
   {

      rnum--;
      twook = 0;

   }

   if (oldrnum > 1)
   {
      
      if (rp[1] < t1 - roottol || rp[1] > t2 + roottol)
	 rnum--;
      else if (twook == 0)
	 rp[0] = rp[1];
   }
   return (rnum);

}

int rotax(const sline &s0,
	  const sline &s1,
	  const sline &s2,
	  const sline &s3,
	  sline *ressl,
	  double curvradmax,
	  double mindist,
	  double planeps)
{
  double n0[3],
         p0[3],
         n1[3],
         p1[3];

  double n[3],
         p[3],
         pi1[3],
         p0i[3],
         p10[3];

  double roots[2];
  double t[2];
  double q0[3],
         q1[3],
         q[3];
  double *rq;
  double u[2][4];
  double v[3];
  double w[2];
  double a,b,c,d0,d1,jacob;
  double mindist2=mindist*mindist;
  int     i,g,h,rno;
  int     i0,i1,i2,i3;
  const sline *ss[4];

  i3 = 0;
  if ((!s0.finite()) || (!s1.finite()) || (!s2.finite()) || (!s3.finite()))
    return -4;   // Now we do not compute infinite intersections



  ss[0] = &s0;
  ss[1] = &s1;
  ss[2] = &s2;
  ss[3] = &s3;

  a = 0.0; 
  i0 = 0; i1 = 1;
  for (i=0; i < 4; ++i) {
    for (g=i+1; g < 4; ++g) {
      b = fabs((*ss[i] % *ss[g]).Dual());
      // Length of transversal * sin angle
      if (b > a) {
	a = b;
	i0 = i;
	i1 = g;
      }
    }
  }

  // Fill out the rest

  i2 = -1;
  for (i=0; i < 4; ++i) {
    if ((i != i0) && (i != i1)) {
      if (i2 < 0)
	i2 = i;
      else
	i3 = i;
    }
  }

  // Compute the coefficients

  ss[i0]->direction(n0);
  ss[i0]->nearp(p0);
  ss[i1]->direction(n1);
  ss[i1]->nearp(p1);


  sub(p10, p1, p0);

  ss[i2]->direction(n);
  ss[i2]->nearp(p);

  sub(pi1, p, p1);
  sub(p0i, p0, p);

  u[0][0]  = triple(p0i,p10,n);     // a   !!! p10 and not pi1 !!!
  u[0][1]  = triple(pi1,n0,n);      // a0
  u[0][2]  = triple(p0i,n1,n);      // a1
  u[0][3]  = triple(n0,n1,n);       // a01

  ss[i3]->direction(n);
  ss[i3]->nearp(p);

  sub(pi1, p, p1);
  sub(p0i, p0, p);

  u[1][0]  = triple(p0i,p10,n);     // b   !!! p10 and not pi1 !!!
  u[1][1]  = triple(pi1,n0,n);
  u[1][2]  = triple(p0i,n1,n);
  u[1][3]  = triple(n0,n1,n);


// std::cout << "a01 = " << u[0][3]  << ",  a0 =  " << u[0][1] << ",  a1 = " << u[0][2] << ", a = " << u[0][0] << "\n";
// std::cout << "b01 = " << u[1][3]  << ",  b0 =  " << u[1][1] << ",  b1 = " << u[1][2] << ", b = " << u[1][0] << "\n";




#ifdef DEBUG_TRACE
  if (debug_axes != 0) {
    dual dd(*ss[i0] % *ss[i1]);
    std::cout << "###i0-i1:" << i0 << ',' << i1 << ':' <<
      dd << "= cos" <<  acos(dd) << std::endl;
  }
#endif
  

  if ((fabs(u[0][3]) < planeps) && (fabs(u[1][3]) < planeps)) {
    // All the input line directions lie in a plane
    // In general we have an infinite line and a finite line solution

    d0 = u[0][2]*u[1][1] - u[1][2]*u[0][1];
    i = 0;

// std::cout <<  "d0 = " <<  d0 << "\n";

    if (fabs(d0) > mindist2) {
      t[0] = (u[0][0]*u[1][2] - u[1][0]*u[0][2])/d0;
      t[1] = (u[1][0]*u[0][1] - u[0][0]*u[1][1])/d0;

// std::cout <<  "t0 = " <<  t[0] <<  "t1 = " <<  t[1] << "\n";

      if ((fabs(t[0]) < curvradmax) && (fabs(t[1]) < curvradmax)) {

	mult(q0, t[0], n0);
	add(q0, q0, p0);

	mult(q1, t[1], n1);
	add(q1, q1, p1);

	ressl[i] = sline(q0, q1, 1);

#ifdef DEBUG_TRACE
  if (debug_axes != 0) {
    std::cout << "*1* ressl[" << i << "]:" << ressl[i] << std::endl;
    std::cout << "*1* i2 dis:" << *ss[i2] % ressl[i] <<
      "    i3 dis:" << *ss[i3] % ressl[i] << std::endl;
    
    for (h=0; h < 2; ++h)
      w[h] = u[h][0] + u[h][1]*t[0] + u[h][2]*t[1] +
	u[h][3]*t[0]*t[1];
    std::cout << "*1* w[0]=" << w[0] <<
      "    w[1]=" << w[1] << std::endl;
    jacob =  u[0][1]*u[1][2] - u[0][2]*u[1][1];
    std::cout << "### Rotax Jacobian:" << jacob << std::endl;
  }
#endif
	i++;
      }
    } // else all straight lines in the plane are solutions

    v_mult(q, n0, n1);

    if (vabs(q) < planeps) {

      // I guess all the (input) straight lines  are parallel
      // This is basically a plane!
      // I think if we are here then it must have not gone in the previous
      // if (fabs(d0) > mindist2). If still, then we shall neglect the fact
      // and return two (!) "infinite" axes which are orthogonal to each other
      // and to the straight lines, thus lying in the plane.
      // This will correspond to a "twofold" translational surface like
      // the plane is.

      makeCoordinateSystem(n0, q0, q1);
      ressl[0] = sline(NULL, q0);
      ressl[1] = sline(NULL, q1);

      // Note, that it is not quite correct since in principle every
      // direction orthognal to 'n0'  should have been assigned to and it is
      // questionable that for other 4-tuple normals of the same plane surface
      // what shall we get. But since the normal 'n0' will be the same there
      // as well and 'makeCoordinateSystem' works in the same way there as
      // well, hopefully we shall get back the same orthogonals.
      // Anyway if something should be signaled extra one can return
      // rno=-2 instead.
      // I still think that the "ordinary" practice will work here as well.
      // Thus:

      rno = 2;

    } else {
      // This is the axis of revolution in the infinity, i.e.
      // a translational direction
      // Note that the 'curvradmax' limit as the distance from the axis
      // is obviously not applied here
      ressl[i] = sline(NULL, q);
      i++;
      rno = i;
    }
  } else {



    for (i=0; i < 3; ++i) {
      v[i] = u[1][3]*u[0][i]  - u[0][3]*u[1][i];
    }

    g = (fabs(v[1]) < fabs(v[2])) ?
      // t1 should be expressed through t0
      // t1 = -v[1]*t0/v[2] - v[0]/v[2]
      2
      :
      1;

    if (fabs(v[g]) < DBL_EPSILON) {
      // Disaster, the only case we can handle is if both u[0][3] and u[1][3]
      // are small. Actually we did this above.
      // parampapam
      if (fabs(v[0]) < DBL_EPSILON) {
	rno = -3;  // Actually I do not know when does it occur.
	           // Infinite number of straight lines.
      } else {
	rno = 0;   // Neither this one. No solution.
      }
    } else {
      d0 = -v[3-g]/v[g];
      d1 = -v[0]/v[g];

      h =  (fabs(u[0][3]) > fabs(u[1][3])) ? 0 : 1;

      a = u[h][3]*d0;
      b = u[h][3]*d1 + u[h][3-g] + u[h][g]*d0;
      c = u[h][g]*d1 + u[h][0];


      rno = v_deg2eqn(a,b,c, -curvradmax, curvradmax, DBL_EPSILON, roots);

      if (rno == -1) {
	q[0] = q[1] = q[2] = 0.0;
	v_mult(q0, n0, n1);
	if (vabs(q0) < DBL_EPSILON) {
	  ss[i2]->direction(n);
	  v_mult(q0, n0, n);
	  ss[i2]->direction(n);
	  v_mult(q1, n1, n);

	  rq = (vabs(q0) < vabs(q1)) ? q1 : q0;

	  if (vabs(rq) < DBL_EPSILON) {
	    rno = -4;
	  } else {
	    ressl[0] = sline(q,rq);
	  }
	} else {
	  ressl[0] = sline(q,q0);
	}
      } else {
	for (i=0; i < rno; ++i ) {
	  t[2-g] = roots[i];
	  t[g-1] = d0*t[2-g] + d1;
	
	  mult(q0, t[0], n0);
	  add(q0, q0, p0);

	  mult(q1, t[1], n1);
	  add(q1, q1, p1);

	  ressl[i] = sline(q0, q1, 1);

#ifdef DEBUG_TRACE
  if (debug_axes != 0) {
    std::cout << "*2* ressl[" << i << "]:" << ressl[i] << std::endl;
    std::cout << "*2* i2 dis:" << *ss[i2] % ressl[i] <<
      "    i3 dis:" << *ss[i3] % ressl[i] << std::endl;
    
    for (h=0; h < 2; ++h)
      w[h] = u[h][0] + u[h][1]*t[0] + u[h][2]*t[1] +
	u[h][3]*t[0]*t[1];
    std::cout << "*2* w[0]=" << w[0] <<
      "    w[1]=" << w[1] << std::endl;
    jacob = v[1]*t[0] - v[2]*t[1] + u[0][1]*u[1][2] - u[0][2]*u[1][1];
    std::cout << "### Rotax Jacobian:" << jacob << std::endl;
  }
#endif
	}
      }
    }
  }


  return rno;

}

void sline::reverse()
{
  if (finite()) {
    dual min1(-1.0,0.0);

    mult(sl, min1, sl);
  } else {
    int  i;

    for (i=0; i < 3; ++i)
      sl[i] = dual(-sl[i].Real(), sl[i].Dual());
  }
}

int guess_axes
(const sline *norms,     // Array of normal lines of a surface
 int          nnorms,    // Number of normals > 4
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
 int         *axsize, // The size info for 'axs' above.
 double       maxsdn,   // Maximum st. deviation in "near" positions
 double     **sdns, // Standard deviation of "near" positions (output)
 int          kslice, // Number of cell slices on a halfsphere
                         // k*(k+1)/2 cells in total
 double       curvradmax, // Maximum values for t0, t1
                                   // only for "finite" axes
 double       mindist, // Minimum meaningful distance
 double       planeps, // Minimum value for angles or such
 int          mincno,  // Min. number of cardinality in a cell
 int        **cnos, // Cardinality of cells (output)
 double       maxsdd,  // Maximum st. deviation in direction
                         // as a proportion to 3.545/k (2sqrt(M_PI)/k)
 double     **sdds) // St. dev of directions (output)
  //
  // Returns the number of axis which have met all the criteria
  // in the order of the cardinality within the cell.
  // Negative returns codes correspond to errors:
  // -1 few input lines
{
  // It will not examine all the \binom^nnorms_4 possibilities, just it
  // makes a linear "shift" around.
  int         error;
  int         i,j,n,m;
  int         nn[2];
  double      dd[2];
  double      dn;
  sline       ressl[2];
  const h_sline *rhsl;
  oh_sline    *res=NULL;
  oh_sline    *x;

  if ((nnorms < 5) || (mincno <= 0)) {
    error = -1;
    goto retret;
  }

  sh_sline::reset_hset(kslice);

  for (i=0; i < nnorms; ++i) {
    n = rotax(norms[i],
	      norms[(i+1) % nnorms],
	      norms[(i+2) % nnorms],
	      norms[(i+3) % nnorms],
	      ressl,
	      curvradmax,
	      mindist,
	      planeps);
    for (j=0; j < n; ++j)
      new sh_sline(ressl[j]);
  }

  m=0;

  for (rhsl = sh_sline::init_cell_sequence();
       rhsl != NULL;
       rhsl = sh_sline::get_next_filled_cell()) {

    rhsl->average(ressl, nn, dd, &dn);
    for (j=0; j < 2; ++j) {
      if ((nn[j] >= mincno) && (dd[j]*kslice < 3.545*maxsdd) &&
	  ((j != 0) || (dn < maxsdn))) {
	x = new oh_sline(ressl[j],
			 nn[j],
			 (j != 0 ? 0.0 : dn),
			 dd[j]);
	x->insert(&res);
	++m;
      }
    }
  }

  if (m > 0) {
    if ((axsize == NULL) || (axsize <= 0)) {
      *axs = new sline[m];
      if (axsize != NULL) {
	*axsize = m;
      }
    } else if (m > *axsize) {
      m = *axsize;              // Output just the first '*axsize' piece.
    }
    if (sdns != NULL)
      *sdns = new double[m];
    if (sdds != NULL)
      *sdds = new double[m];
    if (cnos != NULL)
      *cnos = new int[m];
  } else {
    *axs = NULL;
    if (sdns != NULL)
      *sdns = NULL;
    if (sdds != NULL)
      *sdds = NULL;
    if (cnos != NULL)
      *cnos = NULL;
  }

  x=res;
  for (i=0; i < m; i++) {
    (*axs)[i] = *x;
#ifdef DEBUG_TRACE
  if (debug_axes != 0) {
    std::cout << "$$$ RES line :" << (*axs)[i] <<
      " card:"    << x->get_cn() <<
      " SD(dn):"  << x->get_dn() <<
      " SD(dir):" << x->get_dd() << std::endl;
  }
#endif
    if (sdns != NULL)
      (*sdns)[i] = x->get_dn();
    if (sdds != NULL)
      (*sdds)[i] = x->get_dd();
    if (cnos != NULL)
      (*cnos)[i] = x->get_cn();
    x = (oh_sline *) x->get_next();
  }

  delete res;

  error = m;

retret:
  return error;
}


