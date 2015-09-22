#include "gsbqp.h"

gsBQP::gsBQP(int dim) : elementSupplier(dim)
{ opt = NULL;
   }

gsBQP::~gsBQP()
{ if (opt != NULL) delete opt; 
   }

void gsBQP::init()
{ opt = new solutionBQP(this);
   }

void gsBQP::solve()
{ int i,ok =1,mi;
   double mq,q;

   init();
   mq = opt->quality();

   while (ok)
   { mi = -1;
      for (i = 0; i < n; i++) 
         if (!opt->included(i))  
        { opt->include(i);
           q = opt->quality();
           opt->exclude(i);
           if (q > mq)
           { mq = q;
              mi = i; 
              } 
           }
       if (mi > -1)
       {  opt->include(mi);
           //printf("\nGreedy search: including description %d, quality now %lf",mi,mq);
           } else ok = 0;
       }
   }

int gsBQP::solution(int i)
{ return(opt->included(i)); 
   }

double gsBQP::solution(int* a)
{ int i;
   for (i = 0; i < n; i++) a[i] = opt->included(i);       
   return(opt->quality());
   }


