// model.H

#ifndef LIBSEGMENTOR_MODEL
#define LIBSEGMENTOR_MODEL 1

#include "region.h"
#include "matrix.h"
#include "image.h"

#include <stdio.h>
#include <math.h>

enum { M_NONE, M_XWINDOWS, M_OPENGL };

enum MODELTYPE { CPLANE, CSQ, CSURFACE2, CSPHERE, CCYLINDER, CCONE, CTORUS};

class model
{

public:  
  
  static int image_height, image_width;             // for model drawing
  static double minx,miny,dx,dy;       // for model drawing
  static image *theImage;              // for model drawing
  static int draw_type; // for model drawing;

  // virtual constructors

  virtual model* new_model(region& r) = 0;
  virtual model* improve(region& orig, region& added) = 0;   // calculate new model parameters when
  // region added is added to region orig
  virtual model* improve(model *m, region& orig, region& added) = 0;

  virtual double distance(struct point& p) { return(fabs(signed_distance(p))); } // distance function
  virtual double abs_signed_distance(struct point &p) = 0;
  virtual double signed_distance(struct point &p);
  virtual void draw() = 0;
  virtual void drawGL() = 0;
  virtual void print() = 0;
  virtual void fprint(FILE *f) = 0;
  virtual void save(FILE *f) = 0;                 // save model parameters to file
  virtual int no_param() = 0;                     // number of parameters
  int ih() { return(image_height); }      // image height - necessary for drawing SQ-s
  int iw() { return(image_width); }       // image width
  virtual void get_normal_vector(double &x, double &y, double &z) = 0;
  virtual double model_distance() { return(0.0); }
  virtual MODELTYPE what_model() = 0;
  virtual void rif_write (FILE *f) 
  { if (f) f = f;   // to avoid stupid compiler warnings.
    }
  virtual void parameters(char **name, double *value) = 0;
  virtual void set_parameters(double *value) = 0;
  void describe_model(char *s);
  virtual struct point project(struct point &p) { return(p); }
  virtual void calc_lcs() {}
  virtual vector3D point_to_lcs(struct point& p) 
  { p.x = p.x;  // to prevent stupid compiler warnings
    return(*(new vector3D)); 
    }
  virtual void rotate(symatrix& A) 
  { if (A.el(0,0)) A.el(0,0) = A.el(0,0); // to avoid stupid compiler warnings
    }
  virtual void translate(vector3D& v) 
  { if (v.el(0)) v.el(0) = v.el(0);
    }
  virtual int model_ok() { return(1); }
  };

#endif  




