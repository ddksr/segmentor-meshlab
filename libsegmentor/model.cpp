#include "model.h"
#include "description.h"

int model::image_height = 0;
int model::image_width = 0;
double model::minx = 1e10;
double model::miny = 1e10;
double model::dx = -2e10;
double model::dy = -2e10;     
image *model::theImage = NULL;            
int model::draw_type = M_OPENGL; 

void model::describe_model(char *s)
{ char *str[20];
  double val[20];
  int i,len;

  for (i = 0; i < 20; i++) 
  { str[i] = new char[20];
    str[i][0] = 0;
    }
  parameters((char **)str,val);
  i = len = 0;
  while (str[i][0])
  { sprintf(s+len,"%s=%lf  ",str[i],val[i]);
    len = strlen(s);
    i++;
    }
 
  for (i = 0; i < 20; i++)
    delete [] str[i];

  }

double model::signed_distance(struct point &p)
{ return(abs_signed_distance(p)/description::stdev);
  }


