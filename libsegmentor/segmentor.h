#ifndef SEGMENTOR
#define SEGMENTOR

#define MAX_DESCRIPTION_LISTS 100
#define NUM_OF_MODELS 8

#include "common.h"
#include "image.h"
#include "segmentation.h"


class Segmentor {
 public:
  ~Segmentor();

  static Segmentor* Instance();

  void setUp(RecoverySettings*, image*, Drawer*, ProgressIndicator*, Messaging*);

  void import(MODELTYPE, double*);
  
  Drawer* getDrawer();
  void refreshConfig();
  void placeSeeds();
  void grow();
  void selection();
  void finalSelection();
  void deleteWrong();

  void mergeDescriptions();
  void fillGaps();
  void intersectionRefinement();

  void enableMessaging(bool x) { messagingEnabled = x; }

  segmentation* descriptions;
  int numOfDescriptions;

 private:
  Segmentor();

  image *im, *normals;
  MODELTYPE *modelType;
  
  ShapeSettings* models[NUM_OF_MODELS];
  
  bool initialized;
  bool messagingEnabled;

  int sel, seld; // TODO: what?
  int *dneigh; // TODO: what?
  
  RecoverySettings* conf;
  Messaging* message;
  Drawer* drawer;
  ProgressIndicator* progress;
  static Segmentor *_inst;

  void clear();
  void merge(segmentation &);
  void refine(segmentation &l, region *r = NULL);
  description *bestDescription(region *r, description *old=NULL, double iq = 0.0);
};

#endif



