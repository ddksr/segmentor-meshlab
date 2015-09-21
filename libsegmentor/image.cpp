// image.C

#include "image.h"

#include <math.h>
#include <stdio.h>
#include <stddef.h>

image::image() {}
  
int image::triang_count(int *colors)
{ unsigned int i;
  int k,no = 0,nk;
  for (i = 0; i < number; i++)
  { nk = pix[i].n >> 1;
    for (k = 0; k < nk; k++)
      if (colors[i] == colors[pix[i].neigh[k<<1]] && colors[i] == colors[pix[i].neigh[1+(k<<1)]])
        no++; 
    }             
  return(no/3);
  }  

int image::triang_count()
{ unsigned int i;
  int k,no = 0, nk,nc = 0, si;
  
  for (i = 0; i < number; i++)
  { nk = pix[i].n >> 1;
    nc += nk;
    si = i;
    for (k = 0; k < nk; k++)
      if (si < pix[si].neigh[2*k] && si < pix[si].neigh[2*k+1]) no++;    
    }

  nc /= 3;
  if (nc != no) printf("\nTriangle count difference %d %d",nc,no);
  return(no);

  }
    

void image::bounding_box(double &x1, double &y1, double &z1, double &x2, double &y2, double &z2)
{ x1 = minx;
  y1 = miny;
  z1 = minz;
  x2 = maxx;
  y2 = maxy;
  z2 = maxz;
  }


int operator==(const image& im1, const image& im2)
{ unsigned int i;
  int n1 = im1.points(),n2 = im2.points();

  if (n1 != n2) return(0);

  for (i = 0; i < im1.number; i++)
    if (im1.pixel(i).n != im2.pixel(i).n) return(0);

  return(1);
  }

void image::extremes()
{ unsigned int i;

  for (i = 0; i < number; i++)
  { if (!i || pix[i].x < minx) minx = pix[i].x;
    if (!i || pix[i].x > maxx) maxx = pix[i].x;
    if (!i || pix[i].y < miny) miny = pix[i].y;
    if (!i || pix[i].y > maxy) maxy = pix[i].y;
    if (!i || pix[i].z < minz) minz = pix[i].z;
    if (!i || pix[i].z > maxz) maxz = pix[i].z;
    }
  }

void image::scale(double s)
{ unsigned int i;

  for (i = 0; i < number; i++)
  { pix[i].x *= s;
    pix[i].y *= s;
    pix[i].z *= s;
    }
  minx *= s;
  maxx *= s;
  miny *= s;
  maxy *= s;
  minz *= s;
  maxz *= s;
  }

image::~image()
{ unsigned int i;
  for (i = 0; i < number; i++)
    if (pix[i].n) delete [] pix[i].neigh;
  delete [] pix;
  if (x != NULL) delete [] x;
  if (y != NULL) delete [] y;
  }

double distance(struct point& p1, struct point& p2)
{ return(sqrt((p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y) + (p2.z - p1.z) * (p2.z - p1.z)));
  }

