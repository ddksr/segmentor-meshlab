#include <stddef.h>

#include "segmentor.h"
#include "image.h"

Segmentor* Segmentor::_inst = NULL;

Segmentor* Segmentor::Instance() {
  if (! _inst) {
	_inst = new Segmentor();
  }
  return _inst;
}

void Segmentor::setUp(RecoverySettings* c, image* img) {
  initialized = true;
  conf = c;
  im = img;
  im->calcNormals();
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
