#ifndef LIBSEGMENTOR_COMMON
#define LIBSEGMENTOR_COMMON

#define PI 3.1415926535


#include "model.h"

struct ShapeSettings {
  bool on;
  float dist;
  float err;
};

struct RecoverySettings {
  ShapeSettings plane;
  ShapeSettings superquadric;
  ShapeSettings secondOrderSurface;
  ShapeSettings sphere;
  ShapeSettings cylinder;
  ShapeSettings cone;
  ShapeSettings torus;
  float k2;
  float k3;
  float planarNormalCompatibility;
  int growingSteps;
  int selections;
  int seedSize;
  float regionCompatibility;
  int maxNumOfSeeds;
  float varianceOfNoise;
  float postProcessingMaxError;
  bool isNewDiscrepancy;
  bool usePostProcessing;
};

class Drawer {
 public:
  virtual void prepare(model*) {}
  virtual void draw() {}
  virtual void clear() {}
};

#endif
