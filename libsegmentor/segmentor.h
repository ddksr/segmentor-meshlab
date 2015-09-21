#ifndef SEGMENTOR
#define SEGMENTOR

#define PI 3.1415926535

#include "common.h"
#include "image.h"

class Segmentor {
 public:
  ~Segmentor();
  
  int numOfDescriptors;
  RecoverySettings config;

  static Segmentor* Instance();

  void setUp(RecoverySettings*, image*);

 private:
  Segmentor();

  bool initialized;
  
  image *im;
  
  RecoverySettings* conf;
  static Segmentor *_inst;
  
};

#endif
