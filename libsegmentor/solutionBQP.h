#ifndef LIBSEGMENTOR_SOLUTION_BQP
#define LIBSEGMENTOR_SOLUTION_BQP 1

#include <stdio.h>

class elementSupplier
{ 
   int *elem;
   double *value;

public:
   int n;  

   elementSupplier(int dim);
   virtual ~elementSupplier();

 /* 
** Function to read elemenrs from the matrix of the problem.
 */
   double element(int i, int j);
/*
** Virtual function to read elements from the matrix of the problem.
 * @param i specifies row of the matrix 
 * @param j specifies column of the matrix
 * @return the value of element at position (i,j) in the matrix
 */
   virtual double el_source(int i, int j) = 0;       
  
   };

class solutionBQP
{ elementSupplier *parent;
   int *sol; 
   double *value;
   double q;

public:
   solutionBQP();
   solutionBQP(elementSupplier *p); 
   solutionBQP(const solutionBQP &prev);
   ~solutionBQP();
   solutionBQP& operator=(solutionBQP prev);
   int included(int i) { return(sol[i]); }
   void include(int i);
   void exclude(int i);
   void toggle(int i);
   double quality() { return(q); }
   };

#endif
