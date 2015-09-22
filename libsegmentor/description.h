// description.H

#ifndef LIBSEGMENTOR_DESCRIPTION
#define LIBSEGMENTOR_DESCRIPTION

#include "model.h"
#include "region.h"

#include <stdio.h>

extern int useStatistics;
extern int drawRegions;

class overflowManager
{  int *over;
    int no;
    char filename[200];

public:
   overflowManager(const char *name = "overflow.P", int nmax = 1000);
   virtual ~overflowManager();
  double P(int m, int n);
  int overflowP(int n, double thr);
   
   };

extern overflowManager oM;

class description
{ 
protected:

  int can_grow_var;
  int camera_set;
  double cameraCoordinates[3];

public:
  static double stdev,k2,k3;
  static int msc;
  static int discrepancy;
  model *mmodel;
  region *mregion;
  image* normals;
  region *allowed;

  description() { can_grow_var = 1;  allowed = NULL; }
  description(region& r);
  description(FILE *f,image *im,image *norm, int camera);
  virtual description* new_description(region& r, model *m = NULL) = 0;  // virtual constructor
  virtual void grow(region *ra = NULL);
  double total_error(region& r, int& card); // returns sum of distances + number of points
  double total_error(int& card) { return(total_error(*mregion, card)); }
  
  double error(model *m);  // average sum of distances of points from mregion to model m
  double residual2(region *r = NULL, int usual = 0);         // sum of squared residuals
  double gprobability();
  double MSC(description *d = NULL);
  double MDL(description *d = NULL);
  double BIC(description *d = NULL);
  double AIC(description *d = NULL);
  double error() { return(error(mmodel)); }
  int can_grow() { return(can_grow_var); }
  void set_grow(int growing);
  int no_param() { return(mmodel->no_param()); }
  
  virtual double max_point_distance() = 0;         // maximal allowed distance of point from model
  virtual double max_error() = 0;                  // maximal model error
  virtual void set_max_point_distance(double) = 0;
  virtual void set_max_error(double) = 0;

  virtual int compatible(int i, int j) { return(mmodel->distance(mregion->get_point(i,j)) <= max_point_distance()); } 
  double noise_variance();
  void draw();
  void print();   
  void save(FILE *f); 
  void clear() { can_grow_var = 1; }
  virtual ~description();
  double& cameraX();
  double& cameraY();
  double& cameraZ();    
  int statOK(int *neigh, int *used, int nt, double thr, int &poz, int &neg);
  int check_region(region &r);
  void analyze();
  void setDiscrepancy(int d = 0) { discrepancy = d; }
  static void setMSC(int m = 0) { msc = m; }
  static void setDeviation(double d) { stdev = d; }
  };
  

#endif

