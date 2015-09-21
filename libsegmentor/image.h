// image.H

#ifndef LIBSEGMENTOR_IMAGE
#define LIBSEGMENTOR_IMAGE 1

#include <stddef.h>

enum BOUNDARY {OUTLIER,NO_BOUND,BOUND,IRREG};

struct point
{ double x,y,z,nx,ny,nz;
  int n,valid;
  int *neigh;
  BOUNDARY boundary;
  };

class image
{ 
public:
  struct point *pix;
  unsigned int w,h,number;
  int *x, *y;
  double minx, maxx, miny, maxy, minz, maxz, norm;

  static double camX,camY,camZ;
  double camOrientationX,camOrientationY,camOrientationZ;
  
protected:  
  double max(double c1, double c2) { return(c1 > c2? c1:c2); }
  double min(double c1, double c2) { return(c1 > c2? c2:c1); }
public:
  image();
  
  int inside(unsigned int i) { return (i < number); }
  int inside(unsigned int i, unsigned int j) { return(i < w && j < h); }
  struct point& pixel(unsigned int i, unsigned j) const { return(pix[i+j*w]); }
  struct point& pixel(unsigned int i) const { return(pix[i]); }
  virtual void set_pixel(unsigned i, unsigned j, double xn, double yn, double zn) 
  { int ind = i + j * w;
    pix[ind].x = xn;
    pix[ind].y = yn;
    pix[ind].z = zn; 
    } 
  void set_normal(unsigned i, double xn,double yn, double zn)
  { pix[i].nx = xn;
    pix[i].ny = yn;
    pix[i].nz = zn;
    }
  int neighbours(unsigned int i, unsigned int j = 0) { return(pix[i+j*w].n); }  
  int neighbour(unsigned int i, unsigned int j, unsigned int m) { return(pix[i+j*w].neigh[m]); }
  void add_neighbour(unsigned int i, unsigned int j, unsigned int m)
  { int a = i + j * w;
    if (!(pix[a].n-((pix[a].n >> 3) << 3)))
    { int *temp,k;
      temp = new int[pix[a].n+8];
      for (k = 0; k < pix[a].n; k++) temp[k] = pix[a].neigh[k];
      if (pix[a].n) delete [] pix[a].neigh;
      pix[a].neigh = temp;
      }
    pix[a].neigh[pix[a].n++] = m;  
    }
  
  virtual int width() { return(w); }
  virtual int height() { return(h); }
  int points() const { return(number); }

  virtual int valid_point(struct point& p) = 0;
  
  int triang_count(int *colors);
  int triang_count();
 
  virtual void mark_boundary() {}
  void bounding_box(double &x1, double &y1, double &z1, double &x2, double &y2, double &z2);
  void extremes();
  // friend int operator==(const image& im1, const image& im2);

  void scale(double s);
  virtual ~image();

  virtual image* calcNormals() { return NULL; }
  
  };

  
double distance(struct point& p1,struct point& p2);
  
#endif  




