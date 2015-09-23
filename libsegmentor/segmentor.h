#ifndef SEGMENTOR
#define SEGMENTOR

#include "common.h"
#include "image.h"

class Segmentor {
 public:
  ~Segmentor();
  
  int numOfDescriptors;
  RecoverySettings config;

  static Segmentor* Instance();

  void setUp(RecoverySettings*, image*, Drawer*);
  Drawer* getDrawer();

 private:
  Segmentor();

  bool initialized;
  
  image *im;
  
  RecoverySettings* conf;
  Drawer* drawer;
  static Segmentor *_inst;
  
};

#endif
