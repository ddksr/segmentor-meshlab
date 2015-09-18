#include <stddef.h>

#include "segmentor.h"

Segmentor* Segmentor::_inst = NULL;

Segmentor* Segmentor::Instance() {
  if (! _inst) {
	_inst = new Segmentor();
  }
  return _inst;
}

void Segmentor::setUp(RecoverySettings* c) {
  initialized = true;
  conf = c;
}

Segmentor::Segmentor() {

}

Segmentor::~Segmentor() {

}
