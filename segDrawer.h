#ifndef SEGDRAWER_H
#define SEGDRAWER_H

#include <QList>
#include <common/interfaces.h>

#include "libsegmentor/common.h"

#include "libsegmentor/model.h"
#include "libsegmentor/image.h"
#include "libsegmentor/sq.h"

class MeshlabDrawer : public Drawer {
 public:
  MeshlabDrawer(image*);

  void prepare(model*);
  void clear();
  void draw();

 private:
  image* img;
  QList<model*> models;

  void draw_sq(sq*);

  void draw_test();
};


#endif;
