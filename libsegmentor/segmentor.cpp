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

#include <QDebug>

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

  segmentation::drawer = d;
  segmentation::progress = pi;

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
  
  descriptions = NULL;
  modelType = NULL;
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
  if (!initialized) return;
  drawer->clear();
  refreshConfig();
  
  segmentation ls;
  bool found;
  qDebug() << "Started placing seeds";
  if (conf->seedSize <= 2) return;
  ls.init_list(descriptions[0].segmentationImage, descriptions[0].normals);

  progress->clear(0, NUM_OF_MODELS);
  progress->setProcessName("Placing seeds ... ");
  
  
  for (int i = 0; i < NUM_OF_MODELS; i++) {
	if (! models[i]->on) { continue; }
	for (int j = 0; j < numOfDescriptions; j++) {
	  if (modelType[j] == i) {
		if (!found) {
		  found = true;
		  ls = descriptions[j];
		} else {
		  ls += descriptions[j];
		}
	  }
	}
	descriptions[numOfDescriptions].init_list(im, normals);
	if (found) {
	  qDebug() << "found one";
	  descriptions[numOfDescriptions].place_seeds(conf->seedSize, (MODELTYPE)i, &ls);
	} else {
	  qDebug() << "no previous desc found";
	  descriptions[numOfDescriptions].place_seeds(conf->seedSize, (MODELTYPE)i);
	}
	qDebug() << "placed" << descriptions[numOfDescriptions].n << "seeds";
	if (descriptions[numOfDescriptions].n > 0) { numOfDescriptions++; }

	progress->inc();
  }
  progress->clear();
}

void Segmentor::grow() {
  if (!initialized) return;
  refreshConfig();
  int i,j,c, n = 0, k;

  progress->clear(0, numOfDescriptions * (conf->growingSteps + 1));
  progress->setProcessName("Growing ... ");
  
  for (i = 0; i < numOfDescriptions; i++) n+= descriptions[i].n;
  k = 0;
  for (i = 0; i < conf->growingSteps; i++) {
	for (j = 0; j < numOfDescriptions; j++) {
	  descriptions[j].grow(1);
	}
	progress->inc();
  }

  progress->setProcessName("Checking ... ");

  for (j = 0; j < numOfDescriptions; j++) {
	c = 0;
	for (i = 0; i < descriptions[j].n; i++) {
	  if (descriptions[j].d[descriptions[j].handle[i]]->can_grow()) c++; // TODO: check handle
	}
	progress->inc();
  }

  progress->clear();
  
}

void Segmentor::selection() {
  if (!initialized) return;
  
}

void Segmentor::finalSelection() {
  if (!initialized) return;
}
