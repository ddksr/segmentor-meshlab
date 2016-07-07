#include <math.h>
#include "mathutil.h"

double spow(double x, double y) 
{
  if (x < 0.0)
    return(-pow(-x, y));
  else if (x == 0.0)
    return(0);
  else
    return(pow(x, y));
 
}

