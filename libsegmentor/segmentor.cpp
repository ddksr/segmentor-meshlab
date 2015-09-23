#include <stddef.h>

#include "common.h"
#include "segmentor.h"
#include "image.h"

// TMP
#include "sq.h"
#include "region.h"
// ENDTMP

Segmentor* Segmentor::_inst = NULL;

Segmentor* Segmentor::Instance() {
  if (! _inst) {
	_inst = new Segmentor();
  }
  return _inst;
}

void Segmentor::setUp(RecoverySettings* c, image* img, Drawer* d) {
  initialized = true;
  conf = c;
  im = img;
  drawer = d;

  // TMP
  region r(img);

  for (int j = 0; j < img->height(); j++)
    for (int i = 0; i < img->width(); i++)
      if (img->valid_point(r.get_point(i,j))) r.set_point(i,j);
  
  sq* s = new sq(r);
  s->a1 = 10;
  s->a2 = 10;
  s->a3 = 10;
  s->px = 0;
  s->py = 0;
  s->pz = 0;
  drawer->prepare((model*) s);
  // ENDTMP
  
}

Drawer* Segmentor::getDrawer() {
  if (initialized) {
	return drawer;
  }
  return NULL;
}

Segmentor::Segmentor() {
  initialized = false;
}

Segmentor::~Segmentor() {
  if (initialized) {
	delete im;
	delete conf;
  }
}
