#ifndef SEGMENTOR
#define SEGMENTOR

#define MAX_DESCRIPTION_LISTS 100
#define NUM_OF_MODELS 7

#include "common.h"
#include "image.h"
#include "segmentation.h"


class Segmentor {
 public:
  ~Segmentor();

  static Segmentor* Instance();

  void setUp(RecoverySettings*, image*, Drawer*, ProgressIndicator*, Messaging*);
  
  Drawer* getDrawer();
  void refreshConfig();
  void placeSeeds();
  void grow();
  void selection();
  void finalSelection();

 private:
  Segmentor();

  image *im, *normals;
  MODELTYPE *modelType;
  
  ShapeSettings* models[NUM_OF_MODELS];
	
  segmentation* descriptions;
  int numOfDescriptions;
  bool initialized;

  int sel, seld; // TODO: what?
  int *dneigh; // TODO: what?
  symatrix *A;
  
  RecoverySettings* conf;
  Messaging* message;
  Drawer* drawer;
  ProgressIndicator* progress;
  static Segmentor *_inst;

  void clear();
};

#endif
