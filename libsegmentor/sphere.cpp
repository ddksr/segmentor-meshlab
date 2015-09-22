// sphere.C

#include "sphere.h"
#include "matrix.h"

sphere::sphere(region& r)
{ int i,j,sing,k;
  struct point p,dp[4];
  double x2,y2,z2,pom = 0.0;
  symatrix A(3),B(3);
  vector x(3),y(3);

  for (i = 0; i < 20; i++) data[i] = 0;

  k = 0;
  for (j = r.miny; j <= r.maxy; j++)
    for (i = r.minx; i <= r.maxx; i++)
      if (r.included(i,j))
      { 
        p = r.get_point(i,j);
        dp[k] = p;
        if (k < 3) k++;
         
        x2 = p.x * p.x;
        y2 = p.y * p.y;
        z2 = p.z * p.z;
        data[0] += p.x;
        data[1] += p.y;
        data[2] += p.z;
        data[3] += x2;
        data[4] += p.x * p.y;
        data[5] += p.x * p.z;
        data[6] += y2;
        data[7] += p.y * p.z;
        data[8] += z2;
        data[9] += x2 * p.x;
        data[10] += x2 * p.y;
        data[11] += x2 * p.z;
        data[12] += y2 * p.x;
        data[13] += y2 * p.y;
        data[14] += y2 * p.z;
        data[15] += z2 * p.x;
        data[16] += z2 * p.y;
        data[17] += z2 * p.z;
        data[18]++;
	}
 
  A.el(0,0) = data[0] * data[0] / data[18] - data[3];
  A.el(0,1) = A.el(1,0) = data[0] * data[1] / data[18] - data[4];
  A.el(0,2) = A.el(2,0) = data[0] * data[2] / data[18] - data[5];
  A.el(1,1) = data[1] * data[1] / data[18] - data[6];
  A.el(1,2) = A.el(2,1) = data[1] * data[2] / data[18] - data[7];
  A.el(2,2) = data[2] * data[2] / data[18] - data[8];
  B = A.inverse(sing);
  if (sing)
  { A.el(0,0) = dp[1].x - dp[0].x; A.el(0,1) = dp[1].y - dp[0].y; A.el(0,2) = dp[1].z - dp[0].z;
    A.el(1,0) = dp[2].x - dp[0].x; A.el(1,1) = dp[2].y - dp[0].y; A.el(1,2) = dp[2].z - dp[0].z;
    A.el(2,0) = dp[3].x - dp[0].x; A.el(2,1) = dp[3].y - dp[0].y; A.el(2,2) = dp[3].z - dp[0].z; 
    B = A.inverse(sing);
    x.el(0) = euclidian(dp[1]) - euclidian(dp[0]);
    x.el(1) = euclidian(dp[2]) - euclidian(dp[0]);
    x.el(2) = euclidian(dp[3]) - euclidian(dp[0]);
    if (sing) std::cout << "\nSingularity ecountered while recovering sphere!";
    } else
    { pom = (data[3] + data[6] + data[8]) / data[18];
      x.el(0) = (data[0] * pom - data[9] - data[12] - data[15]) / 2.0;
      x.el(1) = (data[1] * pom - data[10] - data[13] - data[16]) / 2.0;
      x.el(2) = (data[2] * pom - data[11] - data[14] - data[17]) / 2.0;
      }
  y = B * x;
  a = y.el(0);
  b = y.el(1);
  c = y.el(2);
  
  if (sing) radius = sqrt(sqr(a-dp[0].x) + sqr(b - dp[0].y) + sqr(c - dp[0].z));
    else radius = sqrt(a * a + b * b + c * c + pom - 2 * (a * data[0] + b * data[1] + c * data[2]) / data[18]);
  }

sphere::sphere(FILE *f)
{ fread(&a,sizeof(double),1,f);
  fread(&b,sizeof(double),1,f);
  fread(&c,sizeof(double),1,f);
  fread(&radius,sizeof(double),1,f);
  fread(data,sizeof(double),20,f);
  }

model* sphere::improve(region& orig, region& added)
{ sphere *sph = new sphere(*this);
  int i,j,sing,k;
  struct point p,dp[4];
  double x2,y2,z2,pom = 0.0;
  symatrix A(3),B(3);
  vector x(3),y(3);

  k = orig.width();   // to avoid stupid compiler warnings
  k = 0;
  for (j = added.miny; j <= added.maxy; j++)
    for (i = added.minx; i <= added.maxx; i++)
      if (added.included(i,j))
      { 
        p = added.get_point(i,j);
        dp[k] = p;
        if (k < 3) k++;
        x2 = p.x * p.x;
        y2 = p.y * p.y;
        z2 = p.z * p.z;
        sph->data[0] += p.x;
        sph->data[1] += p.y;
        sph->data[2] += p.z;
        sph->data[3] += x2;
        sph->data[4] += p.x * p.y;
        sph->data[5] += p.x * p.z;
        sph->data[6] += y2;
        sph->data[7] += p.y * p.z;
        sph->data[8] += z2;
        sph->data[9] += x2 * p.x;
        sph->data[10] += x2 * p.y;
        sph->data[11] += x2 * p.z;
        sph->data[12] += y2 * p.x;
        sph->data[13] += y2 * p.y;
        sph->data[14] += y2 * p.z;
        sph->data[15] += z2 * p.x;
        sph->data[16] += z2 * p.y;
        sph->data[17] += z2 * p.z;
        sph->data[18]++;
	}
 
  A.el(0,0) = sph->data[0] * sph->data[0] / sph->data[18] - sph->data[3];
  A.el(0,1) = A.el(1,0) = sph->data[0] * sph->data[1] / sph->data[18] - sph->data[4];
  A.el(0,2) = A.el(2,0) = sph->data[0] * sph->data[2] / sph->data[18] - sph->data[5];
  A.el(1,1) = sph->data[1] * sph->data[1] / sph->data[18] - sph->data[6];
  A.el(1,2) = A.el(2,1) = sph->data[1] * sph->data[2] / sph->data[18] - sph->data[7];
  A.el(2,2) = sph->data[2] * sph->data[2] / sph->data[18] - sph->data[8];
  B = A.inverse(sing);
  if (sing)
  { A.el(0,0) = dp[1].x - dp[0].x; A.el(0,1) = dp[1].y - dp[0].y; A.el(0,2) = dp[1].z - dp[0].z;
    A.el(1,0) = dp[2].x - dp[0].x; A.el(1,1) = dp[2].y - dp[0].y; A.el(1,2) = dp[2].z - dp[0].z;
    A.el(2,0) = dp[3].x - dp[0].x; A.el(2,1) = dp[3].y - dp[0].y; A.el(2,2) = dp[3].z - dp[0].z; 
    B = A.inverse(sing);
    x.el(0) = euclidian(dp[1]) - euclidian(dp[0]);
    x.el(1) = euclidian(dp[2]) - euclidian(dp[0]);
    x.el(2) = euclidian(dp[3]) - euclidian(dp[0]);
    } else
    { pom = (sph->data[3] + sph->data[6] + sph->data[8]) / sph->data[18];
      x.el(0) = (sph->data[0] * pom - sph->data[9] - sph->data[12] - sph->data[15]) / 2.0;
      x.el(1) = (sph->data[1] * pom - sph->data[10] - sph->data[13] - sph->data[16]) / 2.0;
      x.el(2) = (sph->data[2] * pom - sph->data[11] - sph->data[14] - sph->data[17]) / 2.0;
      }
  y = B * x;
  sph->a = y.el(0);
  sph->b = y.el(1);
  sph->c = y.el(2);
  
  if (sing) sph->radius = sqrt(sqr(sph->a-dp[0].x) + sqr(sph->b - dp[0].y) + sqr(sph->c - dp[0].z));
    else 
  { sph->radius = sph->a * sph->a + sph->b * sph->b + sph->c * sph->c + pom;
    sph->radius = sph->radius - 2 * (sph->a * sph->data[0] + sph->b * sph->data[1] + sph->c * sph->data[2]) / sph->data[18];
    sph->radius = sqrt(sph->radius);
    }

  return(sph);
  }

void sphere::print()
{ std::cout << "\nSphere: center(" << a << "," << b << "," << c << ") radius:" << radius;
  }

void sphere::fprint(FILE *f)
{ fprintf(f,"4 %lf %lf %lf %lf",a,b,c,radius);
  }

void sphere::save(FILE *f)
{ fwrite(&a,sizeof(double),1,f);
  fwrite(&b,sizeof(double),1,f);
  fwrite(&c,sizeof(double),1,f);
  fwrite(&radius,sizeof(double),1,f);
  fwrite(data,sizeof(double),20,f);
  }

void sphere::get_normal_vector(double &x, double &y, double &z)
{ double dx,dy,dz,norm;
  dx = x - a;
  dy = y - b;
  dz = z - c;
  norm = sqrt(dx*dx+dy*dy+dz*dz);
  x = dx / norm;
  y = dy / norm;
  z = dz / norm;
  }

void sphere::parameters(char **name, double *value)
{ sprintf(name[0],"a"); value[0] = a;
  sprintf(name[1],"b"); value[1] = b;
  sprintf(name[2],"c"); value[2] = c;
  sprintf(name[3],"radius"); value[3] = radius;
  }

void sphere::set_parameters(double *value)
{ a = value[0];
  b = value[1];
  c = value[2];
  radius = value[3];
  }

struct point sphere::project(struct point &p)
{ struct point p1;
  vector3D v,v1;
  double d;

  v1.el(0) = p.x;
  v1.el(1) = p.y;
  v1.el(2) = p.z;
  v.el(0) = p.x - a;
  v.el(1) = p.y - b;
  v.el(2) = p.z - c;
  v = v.normalize();
  d = sqrt((p.x-a)*(p.x-a)+(p.y-b)*(p.y-b)+(p.z-c)*(p.z-c));
  
  v.multiply_col(0,d-radius);
  v1 -= v;

  p1.x = v1.el(0); 
  p1.y = v1.el(1);
  p1.z = v1.el(2);

  return(p1);
  }
