#ifndef SEGDRAWER_H
#define SEGDRAWER_H

#include <QList>
#include <common/interfaces.h>

#include "libsegmentor/common.h"

#include "libsegmentor/model.h"

#include "libsegmentor/sq.h"
#include "libsegmentor/asq.h"
#include "libsegmentor/tsq.h"
#include "libsegmentor/plane.h"
#include "libsegmentor/surface2.h"
#include "libsegmentor/sphere.h"
#include "libsegmentor/torus.h"
#include "libsegmentor/cone.h"
#include "libsegmentor/cylinder.h"

#include "libsegmentor/image.h"

struct Color {
  float r;
  float g;
  float b;
};

class MeshlabDrawer : public Drawer {
 public:
  MeshlabDrawer(image*);

  void prepare(model*, float, float, float);
  void prepare(region*, float, float, float);
  void prepare(model*);
  void prepare(region*);
  
  void clear();
  void draw();

 private:
  image* img;
  QList<model*> models;
  QList<Color*> modelColors;
  QList<region*> regions;
  QList<Color*> regionColors;

  void draw_region(region*, Color*);
  void draw_sq(sq*, Color*);
  void draw_asq(asq*, Color*);
  void draw_plane(plane* m, Color*);
  void draw_surface2(surface2* m, Color*);
  void draw_cone(cone* m, Color*);
  void draw_torus(torus* m, Color*);
  void draw_sphere(sphere* m, Color*);
  void draw_cylinder(cylinder* m, Color*);
  void draw_tsq(tsq*, Color*);
};


#endif;
