// plane.C

#include "plane.h"
#include "matrix.h"

#include <math.h>


extern "C" {
  #include "rcad_ndm_ut.h"
  void djacobi(double **a, int n, double *d, double **v, int *nrot);
  void deigsrt(double *d, double **v, int n);
}

plane::plane(region& r)
{ int i,j;
  struct point p;
  double det;
  double **mdata, **evec,eval[4];
    
  mdata = dmatrix(1,3,1,3);
  evec = dmatrix(1,3,1,3);

  for (i = 0; i < 10; i++) data[i] = 0.0;
        
  for (j = r.miny; j <= r.maxy; j++)
    for (i = r.minx; i <= r.maxx; i++)
      if (r.included(i,j))
      { 
        p = r.get_point(i,j);
        data[0] += p.x;
        data[1] += p.y;
        data[2] += p.z;
        data[3] += p.x * p.x;
        data[4] += p.x * p.y;
        data[5] += p.x * p.z;
        data[6] += p.y * p.y;
        data[7] += p.y * p.z;
        data[8] += p.z * p.z;
        data[9]++;
        }
  
  mdata[1][1] = data[3]-data[0]*data[0]/data[9];
  mdata[2][2] = data[6]-data[1]*data[1]/data[9];
  mdata[3][3] = data[8]-data[2]*data[2]/data[9];
  mdata[1][2] = mdata[2][1] = data[4]-data[0]*data[1]/data[9];
  mdata[1][3] = mdata[3][1] = data[5]-data[0]*data[2]/data[9];
  mdata[2][3] = mdata[3][2] = data[7]-data[1]*data[2]/data[9];
                 
  djacobi(mdata,3,eval,evec,&i);
  deigsrt(eval,evec,3);
  
  det = sqrt(evec[1][3]*evec[1][3]+evec[2][3]*evec[2][3]+evec[3][3]*evec[3][3]);
    
  a = evec[1][3] / det;
  b = evec[2][3] / det;
  c = evec[3][3] / det;

  if (a * r.theImage->camOrientationX + b * r.theImage->camOrientationY + c * r.theImage->camOrientationZ > 0.0)
  { a = -a;
    b = -b;
    c = -c;
    }

  d = (-a*data[0]-b*data[1]-c*data[2])/data[9];

  free_dmatrix(evec,1,3,1,3);     
  free_dmatrix(mdata,1,3,1,3);
  
  }
  
plane::plane(FILE *f)
{ fread(&a,sizeof(double),1,f);
  fread(&b,sizeof(double),1,f);
  fread(&c,sizeof(double),1,f);
  fread(&d,sizeof(double),1,f);  
  fread(data,sizeof(double),10,f);
  }  
  
model* plane::improve(region& orig, region& added)
{ plane *pl = new plane(*this);
  int i,j;
  struct point p;
  double **mdata,**evec,eval[4],det;
  
  mdata = dmatrix(1,3,1,3);
  evec = dmatrix(1,3,1,3);
  

  for (j = added.miny; j <= added.maxy; j++)
    for (i = added.minx; i <= added.maxx; i++)
      if (added.included(i,j))
      { p = added.get_point(i,j);
        pl->data[0] += p.x;
        pl->data[1] += p.y;
        pl->data[2] += p.z;
        pl->data[3] += p.x * p.x;
        pl->data[4] += p.x * p.y;
        pl->data[5] += p.x * p.z;
        pl->data[6] += p.y * p.y;
        pl->data[7] += p.y * p.z;
        pl->data[8] += p.z * p.z;
        pl->data[9]++;
        }
  
  mdata[1][1] = pl->data[3]-pl->data[0]*pl->data[0]/pl->data[9];
  mdata[2][2] = pl->data[6]-pl->data[1]*pl->data[1]/pl->data[9];
  mdata[3][3] = pl->data[8]-pl->data[2]*pl->data[2]/pl->data[9];
  mdata[1][2] = mdata[2][1] = pl->data[4]-pl->data[0]*pl->data[1]/pl->data[9];
  mdata[1][3] = mdata[3][1] = pl->data[5]-pl->data[0]*pl->data[2]/pl->data[9];
  mdata[2][3] = mdata[3][2] = pl->data[7]-pl->data[1]*pl->data[2]/pl->data[9];
                 
  djacobi(mdata,3,eval,evec,&i);
  deigsrt(eval,evec,3);
  
  det = sqrt(evec[1][3]*evec[1][3]+evec[2][3]*evec[2][3]+evec[3][3]*evec[3][3]);
    
  pl -> a = evec[1][3] / det;
  pl -> b = evec[2][3] / det;
  pl -> c = evec[3][3] / det;
  if (pl->a * orig.theImage->camOrientationX + pl->b * orig.theImage->camOrientationY + pl->c * orig.theImage->camOrientationZ > 0.0)
  { pl -> a = -pl -> a;
    pl -> b = -pl -> b;
    pl -> c = -pl -> c;
    }
  pl -> d = (-pl->a*data[0]-pl->b*data[1]-pl->c*data[2])/data[9];

  free_dmatrix(evec,1,3,1,3);     
  free_dmatrix(mdata,1,3,1,3);

  /*  symatrix A(4);
  vector x(4),y(4);
  double det;
  int sing;

  A.el(0,0) = data[3];
  A.el(1,0) = A.el(0,1) = data[4];
  A.el(2,0) = A.el(0,2) = data[5];
  A.el(1,1) = data[6];
  A.el(2,1) = A.el(1,2) = data[7];
  A.el(2,2) = data[8];
  A.el(3,0) = A.el(0,3) = data[0];
  A.el(3,1) = A.el(1,3) = data[1];
  A.el(3,2) = A.el(2,3) = data[2];
  A.el(3,3) = data[9];
  x.el(0) = -data[0]; x.el(1) = -data[1]; x.el(2) = -data[2]; x.el(3) = -data[9];

  y = A.inverse(sing)*x;
  pl->a = y.el(0);
  pl->b = y.el(1);
  pl->c = y.el(2);
  pl->d = y.el(3)+1;
  det = sqrt(pl->a*pl->a+pl->b*pl->b+pl->c*pl->c);  
  pl->a = pl->a/det;
  pl->b = pl->b/det;
  pl->c = pl->c/det;
  pl->d = pl->d/det;
  */

  return(pl);
  }
  
void plane::print()
{ std::cout << "\nPlane: a = " << a << ", b = " << b << ", c = " << c << ", d = " << d;
  }  

void plane::fprint(FILE *f)
{ fprintf(f,"1 %lf %lf %lf %lf",a,b,c,d);
  }
  
void plane::save(FILE *f)
{ fwrite(&a,sizeof(double),1,f);
  fwrite(&b,sizeof(double),1,f);
  fwrite(&c,sizeof(double),1,f);
  fwrite(&d,sizeof(double),1,f);
  fwrite(data,sizeof(double),10,f);
  }  

double plane::model_distance()
{ return(d);
  }
  
void plane::get_normal_vector(double &x, double &y, double &z)
{ x = a;
  y = b;
  z = c;
  }
  
void plane::parameters(char **name, double *value)
{ sprintf(name[0],"a"); value[0] = a;
  sprintf(name[1],"b"); value[1] = b;
  sprintf(name[2],"c"); value[2] = c;
  sprintf(name[3],"d"); value[3] = d;
  }

void plane::set_parameters(double *value)
{ a = value[0];
  b = value[1];
  c = value[2];
  d = value[3];
  }

struct point plane::project(struct point &p)
{ struct point r;
  vector3D v1,v2,v3;

  v1.el(0) = a;
  v1.el(1) = b;
  v1.el(2) = c;
  v2.el(0) = p.x + d * a;
  v2.el(1) = p.y + d * b;
  v2.el(2) = p.z + d * c;
  v3.el(0) = p.x;
  v3.el(1) = p.y;
  v3.el(2) = p.z;
  v1.multiply_col(0,v2.inner(v1));
  v3 = v3 - v1;
  r.x = v3.el(0);
  r.y = v3.el(1);
  r.z = v3.el(2);
  return(r);
  }

void plane::calc_lcs()
{ vector3D v,v1,v2;
  symatrix A(3);

  int i,k,j;
 
  v.el(0) = a;
  v.el(1) = b;
  v.el(2) = c;

  if (fabs(a) < fabs(b))
    if (fabs(a) < fabs(c)) i = 0;
      else i = 2;
    else if (fabs(b) < fabs(c)) i = 1;
      else i = 2;

  v1.el(i) = 0.0;
  k = (1 << ((3-i)&2)) >> 1;
  j = (1 << ((3-i)&1)) >> 1;
  v1.el(k) = -v.el(j);
  v1.el(j) = v.el(k);
  v1 = v1.normalize();
  v2 = v.outer(v1);
  v2 = v2.normalize();
  
  for (i = 0; i < 3; i++)
  { A.el(i,0) = v1.el(i);
    A.el(i,1) = v2.el(i);
    A.el(i,2) = v.el(i);
    }

  LCS = new symatrix(3);
  *LCS = A.inverse(i);
  }

vector3D plane::point_to_lcs(struct point &p)
{ vector3D v,v1;

  v.el(0) = p.x; 
  v.el(1) = p.y;
  v.el(2) = p.z;
  v1 = (*LCS) * v;
  return(v1);
  }

void plane::rotate(symatrix& A)
{ class vector v(4), u(4);
  
  u.el(0) = a;
  u.el(1) = b;
  u.el(2) = c;
  u.el(3) = 1;
  
  v = A * u;

  a = v.el(0);
  b = v.el(1);
  c = v.el(2);

  }

void plane::translate(vector3D& v)
{ vector3D n;
  double dist;

  n.el(0) = a;
  n.el(1) = b;
  n.el(2) = c;

  dist = n.inner(v);
  d -= dist;
  }
