// sq.H

#ifndef LIBSEGMENTOR_SUPERQUADRIC
#define LIBSEGMENTOR_SUPERQUADRIC 1

#include "region.h"
#include "matrix.h"
#include "model.h"

#include <stdio.h>

class sq : public model
{ 
  
  double f(double x, double y, double z) const;
  void sq_draw(int no_hidding = 0);

public:
  double a1,a2,a3,e1,e2,px,py,pz,phi,theta,psi;    // SQ parameters
  hmatrix l_from_g,g_from_l;
  
  
  sq() {}
  sq(region& r);
  sq(FILE *f);
  sq(sq *m, region& r);
  virtual ~sq() { }
  model* new_model(region& r) { return(new sq(r)); }
  model* improve(region& orig, region& added);
  model* improve(model *m, region& orig, region& added);
  double abs_signed_distance(struct point& p);
  void draw() 
  { switch(draw_type)
    { case M_XWINDOWS:
        sq_draw(); 
        break;
      case M_OPENGL:
        drawGL();
        break;
      default:
        break;
      }
    }
  void drawGL();
  void print(); 
  void fprint(FILE *f);
  void save(FILE *f);
  int no_param() { return(11); }
  void get_normal_vector(double &x, double &y, double &z) 
  { if (x > y && y > z) x = x;  // to avoid stupid compiler messages
    }
  
  vector normal(double, double) const;
  vector r(double, double) const;
  vector transform(vector& v);

  MODELTYPE what_model() { return(CSQ); }  

  void parameters(char **name, double *value);
  void set_parameters(double *value);

  double map_eta(double);
  double map_omega(double);

  };


#endif
