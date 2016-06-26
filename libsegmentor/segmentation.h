#ifndef LIBSEGMENTOR_SEGMENTATION
#define LIBSEGMENTOR_SEGMENTATION

// list.H

enum { RSNORMAL, RSGEOMVIEW, RSSTAGES};

#include "description.h"
#include "matrix.h"
#include "image.h"
#include "common.h"

#include "gsbqp.h"
#include <stdio.h>

class segmentation
{ static int dsize;
  static int dn;
  static int usingList;
  double min(double d1, double d2) { return(d1 > d2? d2:d1); }
public:
  static Messaging* message;

  static description **d;
  static double pp_max_err,pp_max_dist;

  static ProgressIndicator *progress;
  static Drawer *drawer;  
  
  int size,n,showtype;
  int *handle;
  image *segmentationImage;
  image *normals;
  double k2,k3;               // normalized to k1
  void throw_away();
  segmentation(int nsize = 4096, int dnsize = 20480);
  segmentation(const segmentation& l);
  segmentation(FILE *f, image *im_1, image *norm_1);
  void load(FILE *f, image *im_1, image *norm_1);
  segmentation& operator=(const segmentation& l);
  segmentation& operator+=(const segmentation& l);
  segmentation operator+(const segmentation l);
  void import(MODELTYPE, double*);
  void place_seeds(int seed_size, MODELTYPE type);
  void place_seeds(int seed_size, MODELTYPE type, segmentation *l);
  void place_seeds(int seed_size, MODELTYPE type, segmentation *l, int maxseeds);
  symatrix& get_sel_matrix();
  void init_list(image *im, image *norm);
  void selection();
  virtual void grow(int iterations);
  void draw(int mode = 0);
  void set_sel_const(double kk2, double kk3);
  static description* create(region& r, MODELTYPE type, model *m = NULL, int camera = 0, double cX = 0.0, double cY = 0.0, double cZ = 0.0);
  friend std::ostream& operator<<(std::ostream& s, segmentation& l);
  void print() { for (int i = 0; i < n; i++) d[handle[i]]->print(); }
  void print_usage() { std::cout << "\nUsing:" << usingList; }
  int ndescriptions() { return(n); }
  void recover_and_select(int seed_size, int fgi, int iter, MODELTYPE type);
  void simpleSelection();
  void simpleRAS(int seed_size, int fgi, int iter, MODELTYPE type);

  void delete_description(int i);    
  void adjacency_graph(int **graph,int limit = 10);
  void exclude_intersections(double ep);
  int join_compatible(int constr);              // joins compatible triangles
  int triangle_free(int i1, int i2, int i3);    // returns 1, if triangle is not own by any region

  void label_points(FILE *f);
  void fprint(FILE *f);
  void model_register(int *k, segmentation &l, symatrix &A);
  void image_transform(symatrix& A);
  void rotate(symatrix &A);
  void translate(vector3D &v);
  void ICP(double dist, segmentation &l,symatrix &A,region &r1, region &r2);
  void ICRP(double dist, segmentation &l,symatrix &A,region &r1, region &r2);
  void ProjectICP(double dist, int nk, int *k1, int *k2, segmentation &l, symatrix &A, region &r2);
  void MedioniICP(double dist, segmentation &l,symatrix &A,region &r1, region &r2);
  void delete_wrong();
  void include(int i) { handle[n++] = i; }
  void add(description *d1) { d[dn] = d1; handle[n++] = dn++; }
  void reset();
  double calculateDeviation(int without = -1);
  double goodness();
  double goodnessOfOne(int k);
  double goodnessOfTwo(description *d1, description *d2, int k = -1, int l = -1);
  double goodnessOfOne(description *d1, int k = -1, int l = -1);
  double intersectionPenalty(int i, int j);
  virtual ~segmentation();  
  };

class gsSimple : public gsBQP
{ segmentation *l;
public:
  gsSimple(segmentation *l1);
  double el_source(int i, int j);
  };


#endif
