#ifndef SEGMENTOR
#define SEGMENTOR

#include "common.h"

class Segmentor {
 public:
  //list *descriptors;

  ~Segmentor();
  
  int numOfDescriptors;
  RecoverySettings config;

  static Segmentor* Instance();

  void setUp(RecoverySettings*);

 private:
  Segmentor();

  bool initialized;
  RecoverySettings* conf;
  static Segmentor *_inst;
  
};

#endif
