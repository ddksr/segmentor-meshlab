#ifndef LIBSEGMENTOR_CHI_MODEL_SELECTION_CRITERIA
#define LIBSEGMENTOR_CHI_MODEL_SELECTION_CRITERIA 1

#include "description.h"
#include "region.h"


class chiMSC
{ double stdev,mean,conf;
   int n,nmg,param,anyway;
   double *mg,*prob;
   double *dist,*undodist;
   int *cls, *undocls, *ndist, *undondist;

   double calcChi();

public:
   chiMSC(double m, double dev, int np, double thr=0.99);
   ~chiMSC();
   chiMSC(const chiMSC &c);

   int classify(int i);
   double absEvaluation(double &thr, int &r);
   double probability();
   double gprobability();
   void clear();
   void include(description *m, region *r, int fixate = 0, int offset = 0);
   void exclude(description *m, region *r, int offset = 0);
   void fixate(description *m, region *r, int offset = 0);
   void setAnyway(int a = 1) { anyway = a; }
   };



#endif
