
#ifndef LIBSEGMENTOR_GREEDY_SEARCH_BQP
#define LIBSEGMENTOR_GREEDY_SEARCH_BQP 1


#include "solutionBQP.h"

#include <Qt> // TODO: migrate
#include <stdio.h>

class gsBQP : public elementSupplier
{ solutionBQP *opt; 

public:
   gsBQP(int dim);
   virtual ~gsBQP();
   
private:
   void init();

public:
   void solve();

   int solution(int i);
   double solution(int* a);
   };


#endif
