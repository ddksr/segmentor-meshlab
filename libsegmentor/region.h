// region.H

#ifndef LIBSEGMENTOR_REGION
#define LIBSEGMENTOR_REGION 1

#include "image.h"
#include "bitvect.h"
 
class region
{ unsigned int w,h;

  int max(int i,int j) { return(i < j ? j:i); }
  int min(int i,int j) { return(i > j ? j:i); }

  
public:
  image *theImage;

  bitvect *b;
  int minx,maxx,miny,maxy;
  
  region(image *im);
  region(unsigned int, unsigned int, image *im);

  region(const region& r);             // copy constructor
  
  region operator=(const region r);
  region& operator|=(const region &r);
  region& operator&=(const region &r);
  region operator~();
  region operator&(const region &r);
  region operator|(const region &r);

  struct point& get_point(unsigned int i, unsigned int j) { return (theImage->pixel(i,j)); }   
  
  int width() const { return(w); }
  int height() const { return(h); }
  int inside(unsigned int i, unsigned int j) { return(i < w && j < h); }
  void adjust(int i, int j)
  { if (i < minx) minx = i;
    if (j < miny) miny = j;
    if (i > maxx) maxx = i;
    if (j > maxy) maxy = j;
    }
  void set_point(int i, int j) 
  { b -> setbit(i,j); 
    adjust(i,j); 
    }
  void set_point(int i) 
  { int k = i / w;
    b -> setbit(i); 
    adjust(i-w*k,k);
    }
  void reset_point(unsigned int i, unsigned int j) { b -> resetbit(i,j); }
  void reset_point(unsigned int i) { b -> resetbit(i); }
  void reset_all();
  void set_all();
  void set_rect(unsigned int , unsigned int , unsigned int , unsigned int);
  int included(unsigned int i, unsigned j) const { return(b -> bit(i,j)); }
  int included(unsigned int i) const { return(b -> bit(i)); }
  
  int point_count() const;
  int triang_count();
  region neighbourhood();
  int neighbourhood(int *neigh, int *used);                // for quick browsing through pixels
  int neighbourhood(int *neigh, int *used, int nt);    // for quick browsing through pixels

  void center_of_mass(double &x, double &y, double &z);  
  int triangle_near(int i1, int i2, int i3);
  
  void print() const;

  double surface() const;
  region boundary();
  void bounding_box(double &x1, double &y1, double &z1, double &x2, double &y2, double &z2);

  int neighbour(region &r);            // returns 1 if r is region's neighbour.
  int connected();                         // returns 1 if region is connected.
  region closest(region &r);          // returns closest subregion of r
  ~region();
  };
  
  
  
#endif







