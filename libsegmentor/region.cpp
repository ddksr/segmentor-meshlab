// region.C

#include "region.h"
#include "image.h"
#include "matrix.h"

#include <iostream>
#include <stdlib.h>

region::region(image *im)
{ if (im != NULL)
  { theImage = im;
    w = theImage -> width();
    h = theImage -> height();
    b = new bitvect(w,h);
    maxx = maxy = -1;
    minx = w + 1;
    miny = h + 1; 
    } else
  { std::cout << "\nA call to region::region() with no image specified!";
    exit(1);
    }
  }
  
region::region(unsigned int wr, unsigned int hr, image *im)
{ theImage = im;
  w = wr;
  h = hr;
  b = new bitvect(w,h);
  maxx = maxy = -1;
  minx = w + 1;
  miny = h + 1;
  }
  
region::region(const region& r)
{ w = r.width();
  h = r.height();
  b = new bitvect(*r.b);
  minx = r.minx;
  maxx = r.maxx;
  miny = r.miny;
  maxy = r.maxy;
  theImage = r.theImage;
  }  
  
region region::operator=(const region r)
{ if (this != &r)
  { delete b;
    w = r.width();
    h = r.height();
    b = new bitvect(*r.b);
    minx = r.minx;
    maxx = r.maxx;
    miny = r.miny;
    maxy = r.maxy;
    theImage = r.theImage;
    }  
  return(*this);
  }
  
region& region::operator|=(const region &r)
{ minx = min(minx,r.minx);
  miny = min(miny,r.miny);
  maxx = max(maxx,r.maxx);
  maxy = max(maxy,r.maxy);
  *b |= *r.b;
  if (theImage != r.theImage) 
    std::cout << "\nWarning: Not the same image in operator |= !";
  return(*this);
  }  
  
region& region::operator&=(const region &r)
{ minx = max(minx,r.minx);
  miny = max(miny,r.miny);
  maxx = min(maxx,r.maxx);
  maxy = min(maxy,r.maxy);
  *b &= *r.b;
  if (theImage != r.theImage)
    std::cout << "\nWarning: not the same image in operator &= !";
  return(*this);
  }  

region region::operator&(const region &r)
{ region a = *this;
  a &= r;
/*  unsigned int i,j;
  int k,f = 0;
  for (i = a.minx; i <= a.maxx; i++)
     for (j = a.miny; j <= a.maxy; j++)
       if (a.included(i,j))
         if (!included(i,j) || !r.included(i,j)) f++;
       
  printf("\n\n\nRegion operator & sucks %d times\n\n\n",f);
  print();
  r.print();
  a.print(); */
  return(a);
  }  

region region::operator|(const region &r)
{ region a = *this;
  a |= r;
  return(a);
  }  
    
region region::operator~()
{ region r = *this;
  r.minx = r.miny = 0;
  r.maxx = r.width() - 1;
  r.maxy = r.height() - 1;
  bitvect c(*b); 
  *r.b = ~c;
  return(r);
  }  
  
void region::set_rect(unsigned int i, unsigned int j, unsigned int di, unsigned int dj)
{ int k,l,mi,mj;
  mi = i + di; mj = j + dj;
  for (k = i; k < mi; k++)
    for (l = j; l < mj; l++) set_point(k,l);
  }
  
int region::point_count() const
{ int i,j,sum = 0;
  
  for (j = miny; j <= maxy; j++)
    for (i = minx; i <= maxx; i++)
      if (included(i,j)) sum++;
      
  return(sum);
  }  

int region::triang_count()
{ int i,j,nt = 0,nn,k,s1,s2,s;
  for (j = miny; j <= maxy; j++)
    for (i = minx; i <= maxx; i++)
      if (included(i,j))
      { s = i + j * width();
        nn = theImage -> neighbours(i,j) >> 1;  
        for (k = 0; k < nn; k++)
        { s1 = theImage -> neighbour(i,j,k<<1);
          s2 = theImage -> neighbour(i,j,1+(k<<1));
          if (included(s1)&&included(s2)&&s<s1&&s<s2) nt++;
          }
        }
        
  return(nt);      
  }  
    
region region::neighbourhood()
{ region r(width(),height(),theImage);
  int i,j,k;
 
  for (j = miny; j <= maxy; j++)
    for (i = minx; i <= maxx; i++)
      if (included(i,j))
        for (k = 0; k < theImage->pixel(i,j).n; k++)
          if (!included(theImage->pixel(i,j).neigh[k])) 
            r.set_point(theImage->pixel(i,j).neigh[k]);
  return(r); 
  }    
  
int region::neighbourhood(int *neigh, int *used)
{ int i,j,k,l = 0,px;
 
   for (j = miny; j <= maxy; j++)
     for (i = minx; i <= maxx; i++)
       if (included(i,j))
         for (k = 0; k < theImage->pixel(i,j).n; k++)
           if (!included(theImage->pixel(i,j).neigh[k])) 
           { px = theImage->pixel(i,j).neigh[k];
              neigh[l++] = px;
              used[px] = 1;
              }
  return(l); 
  }    

int region::neighbourhood(int *neigh, int *used, int nt)
{ int i,k,l = nt,px;

     for (i = 0; i < nt; i++)
         for (k = 0; k < theImage->pixel(neigh[i]).n; k++)
          { px = theImage->pixel(neigh[i]).neigh[k];
             if (!used[px])
             { neigh[l++] = px;
                used[px] = 1;
                }
             }
  return(l); 
  }    

    
void region::center_of_mass(double &x, double &y, double &z)
{ int i,j,n;
  struct point p;
  
  x = y = z = 0;
  n = 0;
  for (j = miny; j <= maxy; j++)
    for (i = minx; i <= maxx; i++)
      if (included(i,j))
      { p = theImage -> pixel(i,j);
        x += p.x;
        y += p.y;
        z += p.z;
        n++;
        }
  if (n)
  { x /= n;
    y /= n;
    z /= n;
    } else
  { x = y = z = 0.0;
    }      
  } 

int region::triangle_near(int i1, int i2, int i3)
{ int k = 0;

  if (included(i1)) k++;
  if (included(i2)) k++;
  if (included(i3)) k++;

  if (k == 1 || k == 2) return(1);
    else return(0);
  }
  
void region::reset_all()
{ 
/*  for (j = miny; j <= maxy; j++)
    for (i = minx; i <= maxx; i++)
      reset_point(i,j); */
  b->reset();
  maxx = maxy = -1;
  miny = theImage->height()+1;
  minx = theImage->width()+1;
  }  
  
void region::set_all()
{ 
/*  for (j = 0; j < height(); j++)
    for (i = 0; i < width(); i++)
      set_point(i,j); */
  b->set();
  maxx = width()-1;
  maxy = height()-1;
  miny = minx = 0;
  }  
  
void region::print() const
{ std::cout << "\nRegion: " << minx << ',' << maxx << ',' << miny << ',' << maxy << " no. of points: "
       << point_count() << " surface: " << surface();
  }  

double region::surface() const
{ int i,j,k;
  struct point p,p2,p3; 
  int kn,i2,i3,i1;
  double d = 0.0;
  vector3D x,y,z;

  for (j = miny; j <= maxy; j++)
    for (i = minx; i <= maxx; i++)
      if (included(i,j))
      { p = theImage->pixel(i,j);
        kn = p.n >> 1;
        i1 = i + theImage->width() * j;  
        for (k = 0; k < kn; k++)
        { i2 = p.neigh[2*k];
          i3 = p.neigh[2*k+1];
          if (i1 < i2 && i1 < i3 && included(i2) && included(i3))
	  { p2 = theImage->pixel(i2);
            p3 = theImage->pixel(i3);
            x.el(0) = p2.x - p.x;
            x.el(1) = p2.y - p.y;
            x.el(2) = p2.z - p.z;
            y.el(0) = p3.x - p.x;
            y.el(1) = p3.y - p.y;
            y.el(2) = p3.z - p.z;
            z = x.outer(y);
            d += z.norm() / 2.0;
	    }
	  } 
        }   
  return(d);
  }

region region::boundary()
{ region r(theImage);
  int i,j,k;
  struct point p;  

  r.reset_all();

  for (j = miny; j <= maxy; j++)
    for (i = minx; i <= maxx; i++)
      if (included(i,j))
      { p = theImage -> pixel(i,j);
        if (p.boundary == BOUND) r.set_point(i,j);
	  else
	{ for (k = 0; k < p.n; k++)
            if (!included(p.neigh[k])) r.set_point(i,j);
	  }
        }

  return (r);
  }

void region::bounding_box(double &x1, double &y1, double &z1, double &x2, double &y2, double &z2)
{ int i,j;
  struct point p;

  theImage->bounding_box(x2,y2,z2,x1,y1,z1);
  for (j = miny; j <= maxy; j++)
    for (i = minx; i <= maxx; i++)
      if (included(i,j))
      { p = theImage -> pixel(i,j);
        if (p.x < x1) x1 = p.x;
        if (p.x > x2) x2 = p.x;
        if (p.y < y1) y1 = p.y;
        if (p.y > y2) y2 = p.y;
        if (p.z < z1) z1 = p.z;
        if (p.z > z2) z2 = p.z;
	}
  }  

int region::neighbour(region &r)
{ int i,j,k,n;

  for (i = minx; i <= maxx; i++)
    for (j = miny; j <= maxy; j++)
       if (included(i,j))
       { n = theImage->neighbours(i,j);
         for (k = 0; k < n; k++)
            if (r.included(theImage->neighbour(i,0,k))) return(1);
         }
  return(0);
  }
  
int region::connected()
{ region r = *this, rn = *this;
  int i,j,bi,bj;
  
  bi = bj = -1;
  for (i = minx; i <= maxx; i++)
    for (j = miny; j <= maxy; j++)
       if (included(i,j))
       { bi = i;
         bj = j;
         }
  r.reset_all();  
  rn.reset_all();
  if (bi >= 0 && bj >= 0)
  { r.set_point(bi,bj);
    rn = r.neighbourhood() & (*this);
    while (rn.point_count())
    { r |= rn;
      rn = r.neighbourhood() & (*this);
      }
    if (r.point_count() == point_count()) return(1);
      else return(0);
    } else return(0);
  }

region region::closest(region &r)
{ region r1 = *this,pr(theImage), nh(theImage);
  int i,j,ok;

  pr.reset_all();
  if (r.point_count())
  { pr = r1 & r;
    while (!pr.point_count())
    { nh = r1.neighbourhood();
      if (!nh.point_count())
      { for (j = r.miny, ok = 0; j <= r.maxy && !ok; j++)
          for (i = r.minx; i <= r.maxx && !ok; i++)
             if (r.included(i,j)) 
             { pr.set_point(i,j);
               ok = 1;
               }
        return(pr);
        }
      r1 |= nh;
      pr = r1 & r;
      }
    }
  return(pr);
  }

region::~region()
{ delete b;
  }  

bool region::is_selected(int i, int j) {
  return is_selected(i);
}

bool region::is_selected(int i) {
  if (theImage->numOfSelectedPoints == 0) return true;
  for (int k = 0; k < theImage->numOfSelectedPoints; k++) {
	if (theImage->selectedPoints[k] == i) return true;
  }
  return false;
}

