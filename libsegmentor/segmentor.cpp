#include <stddef.h>

#include "common.h"
#include "segmentor.h"
#include "image.h"

#include "sphere_d.h"
#include "sq_d.h"
#include "torus_d.h"
#include "planar_d.h"
#include "cone_d.h"
#include "cylinder_d.h"
#include "surface2_d.h"

Segmentor* Segmentor::_inst = NULL;

Segmentor* Segmentor::Instance() {
  if (! _inst) {
	_inst = new Segmentor();
  }
  return _inst;
}

void Segmentor::setUp(RecoverySettings* c, image* img, Drawer* d, ProgressIndicator* pi) {
  clear();
  initialized = true;
  conf = c;
  drawer = d;
  progress = pi;

  im = img;
  normals = img->calcNormals();
  model::theImage = im;
  descriptions = new segmentation[MAX_DESCRIPTION_LISTS];
  modelType = new MODELTYPE[MAX_DESCRIPTION_LISTS];
  for (int i = 0; i < MAX_DESCRIPTION_LISTS; i++) descriptions[i].init_list(im, normals);
  numOfDescriptions = 0;

  refreshConfig();
  models[0] = &conf->plane;
  models[1] = &conf->superquadric;
  models[2] = &conf->secondOrderSurface;
  models[3] = &conf->sphere;
  models[4] = &conf->cylinder;
  models[5] = &conf->cone;
  models[6] = &conf->torus;
}

Drawer* Segmentor::getDrawer() {
  if (initialized) {
	return drawer;
  }
  return NULL;
}

Segmentor::Segmentor() {
  initialized = false;

  im = normals = NULL;
  numOfDescriptions = 0;
  
  dneigh = NULL;
  A = new symatrix(3);
  A->identity();
}

Segmentor::~Segmentor() {
  clear();
  if (initialized) {
	delete conf;
  }
}

void Segmentor::clear() {
  if (im != NULL) {
	delete im;
    im = NULL;
  }
  if (normals != NULL) {
	delete normals;
    normals = NULL;
  }
  if (descriptions != NULL) {
	delete [] descriptions;
    descriptions = NULL;
  }
  if (modelType != NULL) {
	delete [] modelType;
    modelType = NULL;
  }
  numOfDescriptions = 0;
}

void Segmentor::refreshConfig() {
  planar_d::m_dist = conf->plane.dist;
  planar_d::m_err = conf->plane.err;
  sq_d::m_dist = conf->superquadric.dist;
  sq_d::m_err = conf->superquadric.err;
  surface2_d::m_dist = conf->secondOrderSurface.dist;
  surface2_d::m_err = conf->secondOrderSurface.err;
  sphere_d::m_dist = conf->sphere.dist;
  sphere_d::m_err = conf->sphere.err;
  cylinder_d::m_dist = conf->cylinder.dist;
  cylinder_d::m_err = conf->cylinder.err;
  cone_d::m_err = conf->cone.dist;
  cone_d::m_dist = conf->cone.err;
  torus_d::m_err = conf->torus.dist;
  torus_d::m_dist = conf->torus.err;
}

void Segmentor::placeSeeds() {
  segmentation ls;
  bool found;

  if (conf->seedSize <= 2) return;
  ls.init_list(descriptions[0].segmentationImage, descriptions[0].normals);
  for (int i = 0; i < NUM_OF_MODELS; i++) {
	if (! models[i]->on) { return; }
	for (int j = 0; j < numOfDescriptions; j++) {
	  if (modelType[j] == i) {
		if (!found) {
		  found = true;
		  ls = descriptions[j];
		} else {
		  ls += descriptions[j];
		}
	  }
	  descriptions[numOfDescriptions].init_list(im, normals);
	  if (found) {
		descriptions[numOfDescriptions].place_seeds(conf->seedSize, (MODELTYPE)i, &ls);
	  } else {
		descriptions[numOfDescriptions].place_seeds(conf->seedSize, (MODELTYPE)i);
	  }
	  if (descriptions[numOfDescriptions].n > 0) { numOfDescriptions++; }
	}
  }
}

void Segmentor::grow() {
  
}

void Segmentor::selection() {
  
}

void Segmentor::finalSelection() {}
