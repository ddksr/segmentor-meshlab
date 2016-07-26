#include <stddef.h>
#include <sstream>

#include "common.h"
#include "state.h"
#include "segmentor.h"
#include "image.h"

#include "sphere_d.h"
#include "sq_d.h"
#include "asq_d.h"
#include "asq.h"
#include "tsq_d.h"
#include "bsq_d.h"
#include "torus_d.h"
#include "planar_d.h"
#include "cone_d.h"
#include "cylinder_d.h"
#include "surface2_d.h"

#include <QDebug>

LookupState* State::lookup = new LookupState;

Segmentor* Segmentor::_inst = NULL;

Segmentor* Segmentor::Instance() {
  if (! _inst) {
	_inst = new Segmentor();
  }
  return _inst;
}

void Segmentor::setUp(RecoverySettings* c, image* img, Drawer* d, ProgressIndicator* pi, Messaging* mess) {
  clear();
  initialized = true;
  conf = c;
  drawer = d;
  progress = pi;
  message = mess;

  segmentation::drawer = d;
  segmentation::progress = pi;
  segmentation::message = mess;

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
  models[7] = &conf->asq;
  models[8] = &conf->tsq;
  models[9] = &conf->bsq;
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
  drawer = NULL;
  progress = NULL;
  descriptions = NULL;
  modelType = NULL;
  dneigh = NULL;
}

Segmentor::~Segmentor() {
  if (initialized) {
	clear();
	delete conf;
  }
}

void Segmentor::import(MODELTYPE mtype, double *params) {
  if (!initialized) return;
  
  drawer->clear();
  refreshConfig();
  
  descriptions[numOfDescriptions].init_list(im, normals);
  descriptions[numOfDescriptions].import(mtype, params);
  
  numOfDescriptions++;
}

void Segmentor::clear() {
  if (drawer != NULL) {
	drawer->clear();
  }
  if (progress != NULL) {
	progress->clear();
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
  cone_d::m_dist = conf->cone.dist;
  cone_d::m_err = conf->cone.err;
  torus_d::m_dist = conf->torus.dist;  
  torus_d::m_err = conf->torus.err;

  asq_d::m_err = conf->asq.err;
  asq_d::m_dist = conf->asq.dist;
  tsq_d::m_err = conf->tsq.err;
  tsq_d::m_dist = conf->tsq.dist;
  bsq_d::m_err = conf->bsq.err;
  bsq_d::m_dist = conf->bsq.dist;

  description::k2 = conf->k2;
  description::k3 = conf->k3;
  description::stdev = conf->varianceOfNoise;
  segmentation::pp_max_err = conf->postProcessingMaxError;
  
  description::discrepancy = conf->isNewDiscrepancy ? 1 : 0;
  useStatistics = conf->useStatistics ? 1 : 0;

  State::lookup->on = conf->useLookup;

  asq::asqType = conf->useAsqKf ? ASQ_TYPE_SINUS : ASQ_TYPE_TAPERING;
  if (conf->useAsqKf) {
	qDebug() << "use sinus";
  } else {
	qDebug() << "use tapering only";
  }
}

description* Segmentor::bestDescription(region *r, description *old, double iq) {
  int i;
  double q,mq = iq;
  description *d, *md;
  d = md = NULL;

  if (r->point_count() < 8) return(NULL);
  for (int i = 0; i < NUM_OF_MODELS; i++) {
    if (! models[i]->on) {
	  d = segmentation::create(*r,(MODELTYPE)i,(old == NULL? NULL:old->mmodel));
      q = d->MSC();
      if (q > mq) {
		if (md != NULL) delete md;
        md = d;
        mq = q;
	  } else delete d;
	}
	return(md);  
  }
}

void Segmentor::merge(segmentation &l) {
  int i,j,bi,bj,ok = 1,k;
  description *d,*md;
  region r(l.segmentationImage);
  segmentation pl;
  double mq,q,bq,a[16];
  int CM = 20479;

  bi = bj = 0;
  mq = l.goodness();
  while (l.n > 1 && ok) {
	ok = 0;
    md = NULL;
    bq = mq;
    for (i = 0; i < l.n - 1; i++)
       for (j = i+1; j < l.n; j++)
         if (l.d[l.handle[i]]->mregion->neighbour(*l.d[l.handle[j]]->mregion))
         { r = (*l.d[l.handle[i]]->mregion) | (*l.d[l.handle[j]]->mregion);
           d = bestDescription(&r);
           if (d != NULL)
           { a[0] = l.d[l.handle[i]]->residual2();
             a[1] = l.d[l.handle[j]]->residual2();
             a[2] = d->residual2();
             a[3] = l.d[l.handle[i]]->MSC();
             a[4] = l.d[l.handle[j]]->MSC();
             a[5] = d->MSC();
             a[6] = l.d[l.handle[i]]->gprobability();
             a[7] = l.d[l.handle[j]]->gprobability();
             a[8] = d->gprobability();
             //qi = l.goodnessOfOne(i);
             //qj = l.goodnessOfOne(j);
             //qij = l.intersectionPenalty(i,j);
             //q = bq - qi - qj + qij + l.goodnessOfOne(d,i,j);
             l.d[CM] = d;
             pl.n = 0;
             for (k = 0; k < l.n; k++)
                if (k != i && k != j) pl.handle[pl.n++] = l.handle[k];
             pl.handle[pl.n++] = CM;
             q = pl.goodness();
             if (q > mq)
             { mq = q;
               bi = i;
               bj = j;
               printf("\n New maxima found by merging descriptions %d and %d  %lf %lf",i,j,q,mq);
               fflush(stdout);
               if (md != NULL) delete md;
               md = d;
               } else delete d;
             }
           }
    if (md != NULL) {
	  ok = 1;
      delete l.d[l.handle[bi]];
      delete l.d[l.handle[bj]];
      printf("\nMerging descriptions %d %d improving quality from %lf to %lf\n\n",bi,bj,mq,bq);
      fflush(stdout);
      l.d[l.handle[bi]] = md;
      l.d[l.handle[bj]] = NULL;
      l.throw_away();
	}
  }
}

void Segmentor::refine(segmentation &l, region *r) {
  region pr(l.segmentationImage),ri(l.segmentationImage),rj(l.segmentationImage),tri(l.segmentationImage),trj(l.segmentationImage);
  region xi(l.segmentationImage),xj(l.segmentationImage),dr(l.segmentationImage);
  int ok,i,j,k,gp,gn = l.n;
  description *d1, *d2;
  double mq,q;

  ok = 1;
  while (ok) {
	ok = 0;
    for (i = 0; i < gn-1; i++)
      for (j = i+1; j < gn; j++) {
		printf("\nRefining descriptions %d and %d",i,j);
         if (r == NULL) pr = (*l.d[l.handle[i]]->mregion) & (*l.d[l.handle[j]]->mregion);
           else pr = (*r);
         if (pr.point_count()) {
		   if (r == NULL && l.d[l.handle[i]]->mregion->point_count() >= 32 &&
			   l.d[l.handle[j]]->mregion->point_count() >= 32) 
             pr |= pr.neighbourhood();
           
           tri = (*l.d[l.handle[i]]->mregion) & ~pr;
           trj = (*l.d[l.handle[j]]->mregion) & ~pr;
           mq = l.d[l.handle[i]]->MSC() + l.d[l.handle[j]]->MSC() + 
                   2 * l.d[l.handle[i]]->MSC(l.d[l.handle[j]]);
           printf("\n\nComparing descriptions %d and %d, initial quality %lf",i,j,mq);
           for (k = 0; k < 2; k++) {
			 if (!k) ri.reset_all();
			 else rj.reset_all();
             if (!k) rj = pr;
			 else ri = pr;

			 do {
			   xi = tri | ri;
               xj = trj | rj;
               printf("\npoints distribution %d - %d ",ri.point_count(), rj.point_count());
               d1 = bestDescription(&xi);
               d2 = bestDescription(&xj);
               if (d1 != NULL && d2 != NULL) {
				 q = d1->MSC() + d2->MSC() + 2 * d1->MSC(d2);
				 printf(" quality %lf",q);
				 if (q > mq) {
				   printf(" -> new maxima found %lf (%lf)",q,mq);
				   mq = q;
				   delete l.d[l.handle[i]];
				   delete l.d[l.handle[j]];
				   l.d[l.handle[i]] = d1;
				   l.d[l.handle[j]] = d2;
				   // ok = 1;
				 }
			   } else printf(" not adequate descriptions");
               if (!k) {
				 if (ri.point_count()) {
				   dr = ri.neighbourhood() & pr;
				   if (!dr.point_count()) dr = ri.closest(rj);
				   ri |= dr;
				 } else ri = tri.closest(pr);
				 rj = pr & ~ri;
			   } else {
				 if (rj.point_count()) {
				   dr = rj.neighbourhood() & pr;
				   if (!dr.point_count()) dr = rj.closest(ri);
				   rj |= dr;
				 } else rj = trj.closest(pr);
				 ri = pr & ~rj;
			   }
			   if (!k) gp = rj.point_count();
			   else gp = ri.point_count();
			 } while (gp);
		   } // for (k = 0; ...)
		 } // if (pr.point_count)
	  } // for (j = i + 1; ...)
  }
  fflush(stdout);
}

void Segmentor::placeSeeds() {
  if (!initialized) return;

  drawer->clear();
  refreshConfig();
  
  segmentation ls;
  bool found;

  if (conf->seedSize <= 2) return;
  ls.init_list(descriptions[0].segmentationImage, descriptions[0].normals);

  progress->clear(0, NUM_OF_MODELS);
  progress->setProcessName("Placing seeds ... ");

  std::ostringstream stream;
  
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
	  descriptions[numOfDescriptions].place_seeds(conf->seedSize, (MODELTYPE)i, &ls);
	} else {
	  descriptions[numOfDescriptions].place_seeds(conf->seedSize, (MODELTYPE)i);
	}

	stream << "Placed " << descriptions[numOfDescriptions].n << " seeds";
	
	if (descriptions[numOfDescriptions].n > 0) { numOfDescriptions++; }

	progress->inc();
  }
  progress->clear();
  if(messagingEnabled) message->info(stream.str().c_str());
}

void Segmentor::grow() {
  if (!initialized) return;
  refreshConfig();
  int i,j,c, n = 0;

  progress->clear(0, numOfDescriptions * (conf->growingSteps + 1));
  progress->setProcessName("Growing ... ");

  drawer->clear();
  
  for (i = 0; i < numOfDescriptions; i++) n+= descriptions[i].n;

  for (i = 0; i < conf->growingSteps; i++) {
	for (j = 0; j < numOfDescriptions; j++) {
	  descriptions[j].grow(1);
	  
	}
	progress->inc();
  }

  progress->setProcessName("Checking descriptions... ");

  for (j = 0; j < numOfDescriptions; j++) {
	c = 0;
	for (i = 0; i < descriptions[j].n; i++) {
	  if (descriptions[j].d[descriptions[j].handle[i]]->can_grow()) c++;
	  drawer->prepare(descriptions[j].d[descriptions[j].handle[i]]->mmodel);
	  drawer->prepare(descriptions[j].d[descriptions[j].handle[i]]->mregion);
	}
	progress->inc();
  }

  progress->clear();
  
}

void Segmentor::selection() {
  int old;
  if (!initialized || !numOfDescriptions) return;
  refreshConfig();
  std::ostringstream stream;

  drawer->clear();
  progress->clear(0, numOfDescriptions);
  progress->setProcessName("Selecting ... ");
  
  for (int j = 0; j < numOfDescriptions; j++) {
	old = descriptions[j].n;
	descriptions[j].set_sel_const(conf->k2, conf->k3);
	descriptions[j].selection();
	stream << descriptions[j].n << " out of " << old << " descriptions selected.\n";

	for (int i = 0; i < descriptions[j].n; i++) {
	  drawer->prepare(descriptions[j].d[descriptions[j].handle[i]]->mmodel);
	  drawer->prepare(descriptions[j].d[descriptions[j].handle[i]]->mregion);
	}
	progress->inc();
  }
  progress->clear();
  if(messagingEnabled) message->info(stream.str().c_str());
}

void Segmentor::finalSelection() {
  if (!initialized || !numOfDescriptions) return;
  segmentation la;
  int previous, i, j;
  std::ostringstream stream;
  
  refreshConfig();
  
  for (i = 0; i < numOfDescriptions; i++) la += descriptions[i];
  la.set_sel_const(conf->k2, conf->k3);
  previous = la.n;
  la.selection();
  stream << la.n << " out of " << previous << " descriptions selected in final selection. ";

  drawer->clear();
  progress->clear(0, numOfDescriptions);
  progress->setProcessName("Final selection ... ");
  for (j = 0; j < numOfDescriptions; j++) {
	previous = descriptions[j].n;
	descriptions[j].grow(0); // throw away handles that point to NULL

	for (i = 0; i < descriptions[j].n; i++) {
	  drawer->prepare(descriptions[j].d[descriptions[j].handle[i]]->mmodel);
	  drawer->prepare(descriptions[j].d[descriptions[j].handle[i]]->mregion);
	}
	progress->inc();
  }
  progress->clear();
  if(messagingEnabled) message->info(stream.str().c_str());
}

void Segmentor::deleteWrong() {
  if (!initialized || !numOfDescriptions) return;
  std::ostringstream stream;
  refreshConfig();

  drawer->clear();
  progress->clear(0, numOfDescriptions);
  progress->setProcessName("Final selection ... ");
  
  for (int i = 0; i < numOfDescriptions; i++) {
	int prev = descriptions[i].n;
	descriptions[i].delete_wrong();
	stream << (prev - descriptions[i].n) << " out of " << prev << " descriptions not OK and deleted.";

	for (int j = 0; j < descriptions[i].n; j++) {
	  drawer->prepare(descriptions[i].d[descriptions[i].handle[j]]->mmodel);
	  drawer->prepare(descriptions[i].d[descriptions[i].handle[j]]->mregion);
	}
	progress->inc();
  }
  progress->clear();
  if(messagingEnabled) message->info(stream.str().c_str());
}

void Segmentor::mergeDescriptions() {
  if (!initialized || !numOfDescriptions) return;

  std::ostringstream stream;
  refreshConfig();

  drawer->clear();
  progress->clear(0, numOfDescriptions + 1);
  progress->setProcessName("Merging descriptions ... ");

  segmentation la;
  la.init_list(descriptions[0].segmentationImage, descriptions[0].normals);
  la = descriptions[0];
  for (int i = 1; i < numOfDescriptions; i++) {
	la += descriptions[i];
	progress->inc();
  }
  la.throw_away();
  merge(la);
  
  progress->inc();
  for (int i = 0; i < numOfDescriptions; i++) {
	descriptions[i].throw_away(); // added this because it didnt work else
	for (int j = 0; j < descriptions[i].n; j++) {
	  drawer->prepare(descriptions[i].d[descriptions[i].handle[j]]->mmodel);
	  drawer->prepare(descriptions[i].d[descriptions[i].handle[j]]->mregion);
	}
  }
  stream << "Merged";
  progress->clear();
  if(messagingEnabled) message->info(stream.str().c_str());
}

void Segmentor::intersectionRefinement() {
  if (!initialized || !numOfDescriptions) return;

  std::ostringstream stream;
  refreshConfig();

  drawer->clear();
  progress->clear(0, numOfDescriptions + 1);
  progress->setProcessName("Refinement ... ");

  segmentation la;
  la.init_list(descriptions[0].segmentationImage, descriptions[0].normals);
  la = descriptions[0];
  for (int i = 1; i < numOfDescriptions; i++) {
	la += descriptions[i];
	progress->inc();
  }
  la.throw_away();
  refine(la);

  for (int i = 0; i < numOfDescriptions; i++) {
	descriptions[i].throw_away(); // added this because it didnt work else
	for (int j = 0; j < descriptions[i].n; j++) {
	  drawer->prepare(descriptions[i].d[descriptions[i].handle[j]]->mmodel);
	  drawer->prepare(descriptions[i].d[descriptions[i].handle[j]]->mregion);
	}
  }

  stream << "Completed";
  progress->clear();
  if(messagingEnabled) message->info(stream.str().c_str());
}

void Segmentor::fillGaps() {
  if (!initialized || !numOfDescriptions) return;

  std::ostringstream stream;
  refreshConfig();

  drawer->clear();
  progress->clear(0, 1);
  progress->setProcessName("Gap filling ... ");

  segmentation la,lc;
  description *d;
  int i,j,k;

  la.init_list(descriptions[0].segmentationImage,descriptions[0].normals);
  lc.init_list(descriptions[0].segmentationImage,descriptions[0].normals);
  la = descriptions[0];
  for (i = 1; i < numOfDescriptions; i++) la += descriptions[i];
  la.throw_away();                                 
  region r(la.segmentationImage),rn(la.segmentationImage),ro(la.segmentationImage);

  for (i = 0; i < la.n; i++) r |= (*la.d[la.handle[i]]->mregion);
  r = ~r;
  while (r.point_count()) {      
	for (i = r.minx; i <= r.maxx; i++)
	  for (j = r.miny; j <= r.maxy; j++)
		if (r.included(i,j)) {
		  ro.reset_all();
		  ro.set_point(i,j);
		  do {
			rn = ro.neighbourhood() & r;
			ro |= rn;
		  } while (rn.point_count());

		  lc.reset();
		  la.throw_away();
		  for (k = 0; k < la.n; k++)
			if (la.d[la.handle[k]]->mregion->neighbour(ro)) {
				lc.include(la.handle[k]);
			}
		  switch(lc.n) {
		  case 1:
			ro |= *lc.d[lc.handle[0]]->mregion;
			d = bestDescription(&ro);
			if (d != NULL) {
			  delete lc.d[lc.handle[0]];
			  lc.d[lc.handle[0]] = d;
			}
			break;    
		  default:
			refine(lc,&ro);
			break;
		  }
		  r &= ~ro;
		  fflush(stdout); 
		}
  }

  for (i = 0; i < numOfDescriptions; i++) {
	descriptions[i].throw_away(); // added this because it didnt work else
	for (j = 0; j < descriptions[i].n; j++) {
	  drawer->prepare(descriptions[i].d[descriptions[i].handle[j]]->mmodel);
	  drawer->prepare(descriptions[i].d[descriptions[i].handle[j]]->mregion);
	}
  }
  
  stream << "Gaps filled";
  progress->inc();
  progress->clear();
  if(messagingEnabled) message->info(stream.str().c_str());
}
