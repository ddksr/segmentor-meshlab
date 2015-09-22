// Projection plane class definitions

#include "pplane.h"
#include "plane.h"

double point_distance(double x1, double y1, double x2, double y2)
{ return((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
  }


pplane::pplane(region &reg, int *ind)
{ plane *pl;
  vector3D *t,n2,*tt;
  symatrix A(3),B(3);
  int i,j, np = reg.point_count(), k;  
  struct point p;
  double d,dist,minx,miny,maxx,maxy,x[9],y[9],ds[9];
  double m11,m20,m02,phi,avx,avy;


  pl = new plane(reg);
  pl->get_normal_vector(n.el(0),n.el(1),n.el(2));
  dist = pl->model_distance();  

  t = new vector3D[np];
  tt = new vector3D[np];

  k = 0;

  for (i = reg.minx; i <= reg.maxx; i++)
    for (j = reg.miny; j <= reg.maxy; j++)
      if (reg.included(i,j))
      { p = reg.get_point(i,j);
        tt[k].el(0) = t[k].el(0) = p.x;
        tt[k].el(1) = t[k].el(1) = p.y;
        tt[k].el(2) = t[k].el(2) = p.z;
        d = n.inner(t[k])+dist;
        n2 = n;
        n2.multiply_col(0,d);
	t[k] -= n2;
        tt[k] -= n2;                      // makes projections to plane
        k++;
	}

  // find vectors that span the plane

  r1 = t[1] - t[0];   // take difference vector of two points on the plane
  r1.normalize();
  r2 = n.outer(r1);   // make the second vector out of vector product between normal and r1

  // find the transformation matrix to r1,r2,n coordinate system

  for (i = 0; i < 3; i++)
  { A.el(i,0) = r1.el(i);
    A.el(i,1) = r2.el(i);
    A.el(i,2) = n.el(i);
    }

  B = A.inverse(i);

  // transform the points

  for (i = 0; i < np; i++)
  { t[i] = B * t[i];          
    // printf("\n%lf %lf %lf",t[i].el(0),t[i].el(1),t[i].el(2));
    }  

  // get the principal axis of the projected set of points
  m11 = m20 = m02 = avx = avy = 0.0;
  
  // compute the center of mass

  for (i = 0; i < np; i++)
  { avx += t[i].el(0);
    avy += t[i].el(1);
    }
  avx /= np;
  avy /= np;
  
  // compute moments of inertia

  for (i = 0; i < np; i++)
  { m11 += (t[i].el(0)-avx)*(t[i].el(1)-avy);
    m20 += (t[i].el(0)-avx)*(t[i].el(0)-avx);
    m02 += (t[i].el(1)-avy)*(t[i].el(1)-avy);
    }

  m11 /= np;
  m20 /= np;
  m02 /= np;

  // compute angle of rotation

  phi = atan(2*m11/(m20-m02));

  // rotate the plane vectors

  vector3D nv1,nv2,mv1,mv2;

  nv1 = mv1 = r1;
  nv1.multiply_col(0,cos(phi));
  mv1.multiply_col(0,-sin(phi));
  nv2 = mv2 = r2;
  nv2.multiply_col(0,sin(phi));
  mv2.multiply_col(0,cos(phi));

  r1 = nv1 + nv2;
  r2 = mv1 + mv2;  

  // compute the transformation matrix again
    
  for (i = 0; i < 3; i++)
  { A.el(i,0) = r1.el(i);
    A.el(i,1) = r2.el(i);
    A.el(i,2) = n.el(i);
    }

  B = A.inverse(i);

  // transform the points

  for (i = 0; i < np; i++)
  { t[i] = B * tt[i];          
    // printf("\n%lf %lf %lf",t[i].el(0),t[i].el(1),t[i].el(2));
    }  

  // find extreme points;

  minx = maxx = t[0].el(0);
  miny = maxy = t[0].el(1); 
  
  for (i = 1; i < np; i++)
  { if (t[i].el(0) < minx) minx = t[i].el(0);
    if (t[i].el(0) > maxx) maxx = t[i].el(0);
    if (t[i].el(1) < miny) miny = t[i].el(1);
    if (t[i].el(1) > maxy) maxy = t[i].el(1);
    }   

  x[0] = x[3] = x[6] = minx; 
  y[0] = y[1] = y[2] = miny;
  x[1] = x[4] = x[7] = (maxx+minx)/2;
  y[3] = y[4] = y[5] = (maxy+miny)/2;
  x[2] = x[5] = x[8] = maxx;
  y[6] = y[7] = y[8] = maxy; 

  for (j = 0; j < 9; j++)
  { ind[j] = 0;
    ds[j] = point_distance(t[0].el(0),t[0].el(1),x[j],y[j]);
    }

  for (i = 1; i < np; i++)
    for (j = 0; j < 9; j++)
    { d = point_distance(t[i].el(0),t[i].el(1),x[j],y[j]);
      if (d < ds[j])
      { ds[j] = d;
        ind[j] = i;
	}
      }

  // printf("\nExtremes:");
  for (i = 0; i < 9; i++)
    { // printf("\n %lf %lf %lf",t[ind[i]].el(0),t[ind[i]].el(1),t[ind[i]].el(2));
    ind[i]++;
    }

  delete [] tt;
  delete [] t;
  delete pl;

  }


// How to use pplane class
/*
void main()
{ triangulated im("c2");
  region r(&im),r1(&im); 
  int i[10];

  r1 = ~r;
  pplane p(r1,i);
  
  }
*/


