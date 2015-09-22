#include "surface2.h"
#include "matrix.h"

#include <math.h>

surface2::surface2(region& r)
{ struct point p;
  int i,j,no;
  double x2,x3,y2,y3;
  int sing;

  vector x(6),y(6);
  symatrix A(6);

  for (no = i = 0; i < 25; i++) data[i] = 0.0;
        
  for (j = r.miny; j <= r.maxy; j++)
    for (i = r.minx; i <= r.maxx; i++)
      if (r.included(i,j))
      { p = r.get_point(i,j);
        x2 = p.x * p.x;
        x3 = x2 * p.x;
        y2 = p.y * p.y;
        y3 = y2 * p.y;
        data[0] += x3*p.x; data[1] += x3*p.y; data[2] += x2*y2; data[3] += p.x*y3;
        data[4] += y3*p.y; data[5] += x3; data[6] += x2*p.y; data[7] += p.x*y2;
        data[8] += y3; data[9] += x2*p.z; data[10] += p.x*p.y*p.z; data[11] += y2*p.z;
        data[12] += x2; data[13] += p.x*p.y; data[14] += y2; data[15] += p.x*p.z;
        data[16] += p.y*p.z; data[17] += p.x; data[18] += p.y; data[19] += p.z;
        no++;
        }

  data[20] = no;
  y.el(0) = data[9]; y.el(1) = data[10]; y.el(2) = data[11];
  y.el(3) = data[15]; y.el(4) = data[16]; y.el(5) = data[19];

  A.el(0,0) = data[0]; A.el(1,0) = A.el(0,1) = data[1];
  A.el(2,0) = A.el(1,1) = A.el(0,2) = data[2];
  A.el(3,0) = A.el(0,3) = data[5];
  A.el(2,1) = A.el(1,2) = data[3];
  A.el(4,0) = A.el(3,1) = A.el(1,3) = A.el(0,4) = data[6];
  A.el(2,2) = data[4];
  A.el(5,0) = A.el(0,5) = data[12];
  A.el(4,1) = A.el(3,2) = A.el(2,3) = A.el(1,4) = data[7];
  A.el(5,1) = A.el(1,5) = data[13];
  A.el(4,2) = A.el(2,4) = data[8];
  A.el(3,3) = data[12];
  A.el(5,2) = A.el(2,5) = data[14];
  A.el(4,3) = A.el(3,4) = data[13];
  A.el(5,3) = A.el(3,5) = data[17];
  A.el(4,4) = data[14]; A.el(5,4) = A.el(4,5) = data[18]; A.el(5,5) = data[20];

  x = A.inverse(sing)*y;

  a = x.el(0);
  b = x.el(1);
  c = x.el(2);
  d = x.el(3);
  e = x.el(4);
  f = x.el(5);

  }

surface2::surface2(FILE *ff)
{ fread(&a,sizeof(double),1,ff);
  fread(&b,sizeof(double),1,ff);
  fread(&c,sizeof(double),1,ff);
  fread(&d,sizeof(double),1,ff);
  fread(&e,sizeof(double),1,ff);
  fread(&f,sizeof(double),1,ff);
  fread(data,sizeof(double),25,ff);
  }

model* surface2::improve(region& orig, region& added)
{ surface2 *s2 = new surface2(*this);
  struct point p;
  int i,j;
  double x2,x3,y2,y3;
  int sing;
  vector x(6),y(6);
  symatrix A(6);
  
  j = orig.width();
  
  for (j = added.miny; j <= added.maxy; j++)
    for (i = added.minx; i <= added.maxx; i++)
      if (added.included(i,j))
      { p = added.get_point(i,j);
        x2 = p.x * p.x;
        x3 = x2 * p.x;
        y2 = p.y * p.y;
        y3 = y2 * p.y;
        s2->data[0] += x3*p.x; s2->data[1] += x3*p.y; s2->data[2] += x2*y2; s2->data[3] += p.x*y3;
        s2->data[4] += y3*p.y; s2->data[5] += x3; s2->data[6] += x2*p.y; s2->data[7] += p.x*y2;
        s2->data[8] += y3; s2->data[9] += x2*p.z; s2->data[10] += p.x*p.y*p.z; s2->data[11] += y2*p.z;
        s2->data[12] += x2; s2->data[13] += p.x*p.y; s2->data[14] += y2; s2->data[15] += p.x*p.z;
        s2->data[16] += p.y*p.z; s2->data[17] += p.x; s2->data[18] += p.y; s2->data[19] += p.z;
        s2->data[20] = s2->data[20]+1;
        }

  y.el(0) = s2->data[9]; y.el(1) = s2->data[10]; y.el(2) = s2->data[11];
  y.el(3) = s2->data[15]; y.el(4) = s2->data[16]; y.el(5) = s2->data[19];

  A.el(0,0) = s2->data[0]; A.el(1,0) = A.el(0,1) = s2->data[1];
  A.el(2,0) = A.el(1,1) = A.el(0,2) = s2->data[2];
  A.el(3,0) = A.el(0,3) = s2->data[5];
  A.el(2,1) = A.el(1,2) = s2->data[3];
  A.el(4,0) = A.el(3,1) = A.el(1,3) = A.el(0,4) = s2->data[6];
  A.el(2,2) = s2->data[4];
  A.el(5,0) = A.el(0,5) = s2->data[12];
  A.el(4,1) = A.el(3,2) = A.el(2,3) = A.el(1,4) = s2->data[7];
  A.el(5,1) = A.el(1,5) = s2->data[13];
  A.el(4,2) = A.el(2,4) = s2->data[8];
  A.el(3,3) = s2->data[12];
  A.el(5,2) = A.el(2,5) = s2->data[14];
  A.el(4,3) = A.el(3,4) = s2->data[13];
  A.el(5,3) = A.el(3,5) = s2->data[17];
  A.el(4,4) = s2->data[14]; A.el(5,4) = A.el(4,5) = s2->data[18]; A.el(5,5) = s2->data[20];

  x = A.inverse(sing)*y;

  s2->a = x.el(0);
  s2->b = x.el(1);
  s2->c = x.el(2);
  s2->d = x.el(3);
  s2->e = x.el(4);
  s2->f = x.el(5);

  return(s2);
  }

void surface2::print()
  { std::cout << "\n2nd order surface: a = " << a << ", b = " << b << ", c = " << c << ", d = " << d 
	 << ",e = " << e << ",f = " << f;
  }  

void surface2::fprint(FILE *ff)
{ fprintf(ff,"3 %lf %lf %lf %lf %lf %lf",a,b,c,d,e,f);
  }
  
void surface2::save(FILE *ff)
{ fwrite(&a,sizeof(double),1,ff);
  fwrite(&b,sizeof(double),1,ff);
  fwrite(&c,sizeof(double),1,ff);
  fwrite(&d,sizeof(double),1,ff);
  fwrite(&e,sizeof(double),1,ff);
  fwrite(&f,sizeof(double),1,ff);
  fwrite(data,sizeof(double),25,ff);
  }  
  
void surface2::parameters(char **name, double *value)
{ sprintf(name[0],"a"); value[0] = a;
  sprintf(name[1],"b"); value[1] = b;
  sprintf(name[2],"c"); value[2] = c;
  sprintf(name[3],"d"); value[3] = d;
  sprintf(name[4],"e"); value[4] = e;
  sprintf(name[5],"f"); value[5] = f;
  }

void surface2::set_parameters(double *value)
{ a = value[0];
  b = value[1];
  c = value[2];
  d = value[3];
  e = value[4];
  f = value[5];
  }
