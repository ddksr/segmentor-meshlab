#include "solutionBQP.h"

elementSupplier::elementSupplier(int dim)
{ int i,m;
   n = dim;
   
   m = n * n;
   elem = new int[m];
   if (elem != NULL)
   { value = new double[m];
      if (value != NULL) 
        for (i = 0; i < m; i++) elem[i] = 0;
           else
           { delete [] elem;
              elem = NULL;
              }
      } else value = NULL;
   }

elementSupplier::~elementSupplier()
{ if (elem != NULL) delete [] elem;
   if (value != NULL) delete [] value;
   }

double elementSupplier::element(int i, int j)
{ int k;
   if (elem == NULL) return(el_source(i,j));
     else
    { k = i + j * n;
       if (!elem[k])
      { elem[k] = 1;
         value[k] = el_source(i,j);
         }
      return(value[k]); 
      }
   }

solutionBQP::solutionBQP()
{ sol = NULL;
   value = NULL;
   }

solutionBQP::solutionBQP(elementSupplier *p)
{ int i;
   parent  = p;
   sol = new int[parent->n];
   value = new double[parent->n];
   q = 0.0;

   for (i = 0; i < parent->n; i++)
   { sol[i] = 0;
      value[i] = parent->element(i,i);
      }
   }

solutionBQP::solutionBQP(const solutionBQP &prev)
{ int i;
   parent = prev.parent;

   sol = new int[parent->n];
   value = new double[parent->n];
 

   q = prev.q;

   for (i = 0; i < parent->n; i++)
   { sol[i] = prev.sol[i];
      value[i] = prev.value[i];
      }
   }

solutionBQP::~solutionBQP()
{ delete [] sol;
   delete [] value;
   }

solutionBQP& solutionBQP::operator=(solutionBQP prev)
{ int i;
   parent = prev.parent;

   if (sol != NULL) delete [] sol;
   if (value != NULL) delete [] value;

   sol = new int[parent->n];
   value = new double[parent->n];

   q = prev.q;

   for (i = 0; i < parent->n; i++)
   { sol[i] = prev.sol[i];
      value[i] = prev.value[i];
      }

  return(*this);
  }

void solutionBQP::include(int i)
{ int j;
   if (!sol[i] && i >= 0 && i < parent->n)
   { for (j = 0; j < parent->n; j++) 
         if (j != i) value[j] += parent->element(i,j) + parent->element(j,i);
      q += value[i];
      sol[i] = 1;
      }
   }

void solutionBQP::exclude(int i)
{ int j;
   if (sol[i] && i >= 0 && i < parent->n)
   { for (j = 0; j < parent->n; j++)
         if (j != i) value[j] -= parent->element(i,j) + parent->element(j,i);
      q -= value[i];
      sol[i] = 0;
      }
  }

void solutionBQP::toggle(int i)
{  if (sol[i]) exclude(i);
      else include(i);
   }

