// cylinder.C

#include "cylinder.h"
#include "matrix.h"

#include <iostream>
#include <math.h>


extern "C" {
  #include "rcad_ndm_ut.h"
  #include "vector_def.h"
  #include "num_anal.h"

  int cyl_fit(double **list_x, double *list_y, double **norm, int no, double *param, double *chisq);
}

#define TRUE 1
#define FALSE 0

#define RIDICULOUS_NUMBER -1000000000


cylinder::cylinder(region &reg) {
 
  #ifndef MAX_LEN
    #define MAX_LEN 20000
  #endif
  
  // interface to C code and Numerical Recipes

  double **list_x = dmatrix(1,MAX_LEN,1,2);
  double *list_y  = dvector(1,MAX_LEN);
  double *param = dvector(1,13);
  double chisq;
 
  double **list_norm = dmatrix(1,MAX_LEN,1,3);


  int x, y, min_x,min_y,max_x,max_y, no;
  struct point point;
 
  no = 0;

  min_x = reg.minx;
  min_y = reg.miny;
  max_x = reg.maxx;
  max_y = reg.maxy;
  
  
       
 
   for(x = min_x; x <= max_x; x++)
   for(y = min_y; y <= max_y; y++)
     if (reg.included(x, y))
      if (no < MAX_LEN) { 
	  no++;
 	  point = reg.get_point(x, y);
 	  list_x[no][1] = point.x;
 	  list_x[no][2] = point.y;
	  list_y[no] = point.z; 

          // point = reg.norm(x,y);

          list_norm[no][1] = point.nx;
 	  list_norm[no][2] = point.ny;
	  list_norm[no][3] = point.nz;

	  
	  
	 
 	}	
 	else
  	  std::cout << "Over the limit of Cylinder points";

  param[1] = RIDICULOUS_NUMBER;



   if( cyl_fit(list_x, list_y, list_norm, no, param, &chisq) == 2)
 
   { /* Gabor rep */
    
     a[1] = param[1];
     a[2] = param[2];
     a[3] = param[3];
     a[4] = param[4];
     a[5] = param[5];
     a[11] = param[11];
     a[12] = param[12];
     a[13] = param[13];
     rep = 2;
     adjust_z(reg);
}
   
else


   if( cyl_fit(list_x, list_y, list_norm, no, param, &chisq) == 1)
 
   { /* not a cylinder maintain a general quadric rep */
     a[1] = param[1];
     a[2] = param[2];
     a[3] = param[3];
     a[4] = param[4];
     a[5] = param[5];
     a[6] = param[6];
     a[7] = param[7];
     a[8] = param[8];
     a[9] = param[9];
     a[10]= param[10];
     a[11] = param[11];
     a[12] = param[12];
     a[13] = param[13];
     rep = 1;
     
     /* force big error */
    } else 
  {
  rep = 0;
  norm.x = param[1];
  norm.y = param[2];
  norm.z = param[3];
  disp.x = param[4];
  disp.y = param[5];
  disp.z = param[6];
  r = param[7];
  
  //  cout << "SPHERE NEW\na = " << a << " b = " << b << " c = " << c << " r = " << r << '\n';
  }



  free_dmatrix(list_norm,1,MAX_LEN,1,3);
  free_dvector(param,1,13);

  free_dvector(list_y,1,MAX_LEN);
  free_dmatrix(list_x,1,MAX_LEN,1,2);
}

cylinder::cylinder(region &reg, cylinder &init) {
 
  #ifndef MAX_LEN
    #define MAX_LEN 50000
  #endif

  // interface to C code and Numerical Recipes

  double **list_x = dmatrix(1,MAX_LEN,1,2);
  double *list_y  = dvector(1,MAX_LEN);
  double *param = dvector(1,13);
  double chisq;

   double **list_norm  = dmatrix(1,MAX_LEN,1,3);


  int x, y, min_x,min_y,max_x,max_y, no;
  struct point point;

  no = 0;
  
  min_x = reg.minx;
  min_y = reg.miny;
  max_x = reg.maxx;
  max_y = reg.maxy;

if (init.rep == 2)
     { param[1] = init.a[1];
     param[2] = init.a[2];
     param[3] = init.a[3];
     param[4] = init.a[4];
     param[5] = init.a[5];

     param[11] = init.a[11];
     param[12] = init.a[12];
     param[13] = init.a[13];

}
else


   if (init.rep == 0)
     { param[1] = init.a[1];
     param[2] = init.a[2];
     param[3] = init.a[3];
     param[4] = init.a[4];
     param[5] = init.a[5];
     param[6] = init.a[6];
     
     param[7] = init.a[7];
     param[8] = init.a[8];
     param[9] = init.a[9];
     param[10] = init.a[10];
      param[11] = init.a[11];
     param[12] = init.a[12];
     param[13] = init.a[13];
    }     
   
   
   else {
   norm.x = init.norm.x;
  norm.y = init.norm.y;
  norm.z = init.norm.z;
  disp.x = init.disp.x;
  disp.y = init.disp.y;
  disp.z = init.disp.z;
  r = init.r;
  
  }
     
 
  for(x = min_x; x <= max_x; x++)
   for(y = min_y; y <= max_y; y++)
      if (reg.included(x, y))
	if (no < MAX_LEN) { 
	  no++;
	  point = reg.get_point(x, y);
	  list_x[no][1] = point.x;
	  list_x[no][2] = point.y;
	  list_y[no] = point.z;

          // point = reg.norm(x,y);

          list_norm[no][1] = point.nx;
 	  list_norm[no][2] = point.ny;
	  list_norm[no][3] = point.nz;

	 
	}	
	else
	  std::cout << "Over the limit of Cylinder points";

          
if( cyl_fit(list_x, list_y, list_norm, no, param, &chisq) == 2)
 
   { 
     a[1] = param[1];
     a[2] = param[2];
     a[3] = param[3];
     a[4] = param[4];
     a[5] = param[5];
     a[11] = param[11];
     a[12] = param[12];
     a[13] = param[13];

     rep = 2;
     adjust_z(reg);
}

else

 if( cyl_fit(list_x, list_y, list_norm, no, param, &chisq) == 1)
 
   { 
     a[1] = param[1];
     a[2] = param[2];
     a[3] = param[3];
     a[4] = param[4];
     a[5] = param[5];
     a[6] = param[6];
     a[7] = param[7];
     a[8] = param[8];
     a[9] = param[9];
     a[10]= param[10];
      a[11] = param[11];
     a[12] = param[12];
     a[13] = param[13];
     rep = 1;
     
    }
 else 
  {
  rep = 0;
  norm.x = param[1];
  norm.y = param[2];
  norm.z = param[3];
  disp.x = param[4];
  disp.y = param[5];
  disp.z = param[6];
  r = param[7];
  
  //  cout << "SPHERE NEW\na = " << a << " b = " << b << " c = " << c << " r = " << r << '\n';
  }

  free_dmatrix(list_norm,1,MAX_LEN,1,3);
  free_dvector(param,1,13);

  free_dvector(list_y,1,MAX_LEN);
  free_dmatrix(list_x,1,MAX_LEN,1,2);
}

/*
cylinder::cylinder(region &reg, cylinder& init) {
 
  #ifndef MAX_LEN
    #define MAX_LEN 20000
  #endif

  // interface to C code and Numerical Recipes

  double **list_x = dmatrix(1,MAX_LEN,1,2);
  double *list_y  = dvector(1,MAX_LEN);
  double *param = dvector(1,13);
  double chisq;
 
  double **list_norm = dmatrix(1,MAX_LEN,1,3);


  int x, y, min_x,min_y,max_x,max_y, no;
  struct point point;

  no = 0;

  min_x = reg.minx;
  min_y = reg.miny;
  max_x = reg.maxx;
  max_y = reg.maxy;
  
  
       
 
   for(x = min_x; x <= max_x; x++)
   for(y = min_y; y <= max_y; y++)
     if (reg.included(x, y))
      if (no < MAX_LEN) { 
	  no++;
 	  point = reg.get_point(x, y);
 	  list_x[no][1] = point.x;
 	  list_x[no][2] = point.y;
	  list_y[no] = point.z; 

          // point = reg.norm(x,y);

          list_norm[no][1] = point.nx;
 	  list_norm[no][2] = point.ny;
	  list_norm[no][3] = point.nz;

	  
	  
	 
 	}	
 	else
  	  std::cout << "Over the limit of Cylinder points";

  param[1] = RIDICULOUS_NUMBER;

   pmem();
   if( cyl_fit(list_x, list_y, list_norm, no, param, &chisq) == 2)
 
     { // Gabor rep
     pmem();
     a[1] = param[1];
     a[2] = param[2];
     a[3] = param[3];
     a[4] = param[4];
     a[5] = param[5];
     a[11] = param[11];
     a[12] = param[12];
     a[13] = param[13];
     rep = 2;
}
   
else


   if( cyl_fit(list_x, list_y, list_norm, no, param, &chisq) == 1)
 
     { // not a cylinder maintain a general quadric rep
     a[1] = param[1];
     a[2] = param[2];
     a[3] = param[3];
     a[4] = param[4];
     a[5] = param[5];
     a[6] = param[6];
     a[7] = param[7];
     a[8] = param[8];
     a[9] = param[9];
     a[10]= param[10];
     a[11] = param[11];
     a[12] = param[12];
     a[13] = param[13];
     rep = 1;
     
     // force big error
    } else 
  {
  rep = 0;
  norm.x = param[1];
  norm.y = param[2];
  norm.z = param[3];
  disp.x = param[4];
  disp.y = param[5];
  disp.z = param[6];
  r = param[7];
  
  //  cout << "SPHERE NEW\na = " << a << " b = " << b << " c = " << c << " r = " << r << '\n';
  }


  free_dmatrix(list_x,1,MAX_LEN,1,2);


  free_dmatrix(list_norm,1,MAX_LEN,1,3);

  free_dvector(list_y,1,MAX_LEN);
  free_dvector(param,1,13);
}
*/
// End addition 18.2.1998


cylinder::cylinder(FILE *f) {
 
  //  fread((char *) &disp.x, 1, sizeof(double), f);
  // fread((char *) &disp.y, 1, sizeof(double), f);
  // fread((char *) &disp.z, 1, sizeof(double), f);
  // fread((char *) &norm.x, 1, sizeof(double), f);
  // fread((char *) &norm.y, 1, sizeof(double), f);
  // fread((char *) &norm.z, 1, sizeof(double), f);
  // fread((char *) &r, 1, sizeof(double), f);
  
  rep = 2;
  fread(a,14,sizeof(double),f);
  fread(&zmin,1,sizeof(double),f);
  fread(&zmax,1,sizeof(double),f); 
}

model* cylinder::improve(region& orig, region& added)
{ region r = orig | added;
  return(new cylinder(r));
  }

model* cylinder::improve(model *m, region& orig, region& added)
{ region r = orig | added;
  return(new cylinder(r,*(cylinder *)m));
  }


/*
double cylinder::distance(struct point& p)
{  double X,Y,Z;
   double r_est, interm;
   
if (rep == 2)
    {   Gabor distance approx eqn 3,2   EOC
        double *pl  = dvector(1,3);
        double *acyl  = dvector(1,3);
        double *r0  = dvector(1,3);
        double *p2 = dvector(1,3);
        double pdota, pdotr0;

        pl[1] = p.x - a[11];
        pl[2] = p.y - a[12];
        pl[3] = p.z - a[13];

        r0[1] = cos(a[2])*sin(a[3]);
        r0[2] = sin(a[2])*sin(a[3]);
        r0[3] = cos(a[3]);

        acyl[1] = cos(a[2])*cos(a[3])*cos(a[4])  - sin(a[2])*sin(a[4]);
        acyl[2] = sin(a[2])*cos(a[3])*cos(a[4])  + cos(a[2])*sin(a[4]);
        acyl[3] = -sin(a[3])*cos(a[4]);


        pdota = dot(pl,acyl);
        p2[1] = p.x;
        p2[2] = p.y;
        p2[3] = p.z;
        pdotr0 = dot(p2,r0);

         r_est = a[5]*(modsq(p2) - 2*a[1]*pdotr0 - pdota*pdota + a[1]*a[1])/2.0 + a[1] - pdotr0;
   
         free_dvector(pl,1,3);
         free_dvector(acyl,1,3);
         free_dvector(r0,1,3);
         free_dvector(p2,1,3);

         return  fabs(r_est);

         
    }


  if (rep == 0)
    {
    interm = sqrt(2.0);
    
    X = p.x - a[11];
    Y = p.y - a[12];
    Z = p.z - a[13];
    
    r_est = a[1]*X*X + a[2]*Y*Y + a[3]*Z*Z + interm*a[4]*X*Y + interm*a[5]*X*Z + interm*a[6]*Y*Z + a[7]*X + a[8]*Y + a[9]*Z + a[10];
    
    return fabs(r_est); 
    
    
    }
  
  
  else {

  X = p.x - disp.x;
  Y = p.y - disp.y;
  Z = p.z - disp.z;
  
  interm = norm.y*(norm.x*Y - norm.y*X) + norm.z*(norm.x*Z -norm.z*X);
  
  r_est = interm*interm;
  
  
  interm = norm.x*(norm.y*X - norm.x*Y) + norm.z*(norm.y*Z - norm.z*Y);
  
  r_est += interm*interm;
  
  interm = norm.x*(norm.z*X - norm.x*Z) + norm.y*(norm.z*Y - norm.y*Z);
  
  r_est += interm*interm;
  
  r_est = sqrt(r_est);  
  
       return(0.0);

    return fabs(r_est - r); 
   }
   
} */

double cylinder::abs_signed_distance(struct point& p)
{ vector3D  axis,n,r,v;
  double d = a[1] + 1 / a[5],result,x;

  axis.el(0) = cos(a[2])*cos(a[3])*cos(a[4]) - sin(a[2])*sin(a[4]);
  axis.el(1) = sin(a[2])*cos(a[3])*cos(a[4]) + cos(a[2])*sin(a[4]);
  axis.el(2) = -sin(a[3])*cos(a[4]);

  // cout << "\nAxis: " << axis;

  n.el(0) = d * cos(a[2])*sin(a[3]);
  n.el(1) = d * sin(a[2])*sin(a[3]);
  n.el(2) = d * cos(a[3]);

  // cout << "\n(rho+1/k)*n: " << n;

  r.el(0) = p.x-a[11];
  r.el(1) = p.y-a[12];
  r.el(2) = p.z-a[13];

  // cout << "\nr: " << r;

  x = r.inner(axis);

  axis.el(0) *= x;
  axis.el(1) *= x;
  axis.el(2) *= x;

  v = n + axis - r;
  
  // cout << "\nv: " << v;

  result = v.norm() - fabs(1/a[5]);
  return(result);
  }

void cylinder::save(FILE *f) { 

  // if ( rep == 1) {

  // fwrite((char *) &disp.x, 1, sizeof(double), f);
  // fwrite((char *) &disp.y, 1, sizeof(double), f);
  // fwrite((char *) &disp.z, 1, sizeof(double), f);
  // fwrite((char *) &norm.x, 1, sizeof(double), f);
  // fwrite((char *) &norm.y, 1, sizeof(double), f);
  // fwrite((char *) &norm.z, 1, sizeof(double), f);
  // fwrite((char *) &r, 1, sizeof(double), f);
  // }

  fwrite(a,14,sizeof(double),f);
  fwrite(&zmin,1,sizeof(double),f); 
  fwrite(&zmax,1,sizeof(double),f);
}

void cylinder::adjust_z(region& reg)
{ vector3D axis,n,pt;
  struct point p;
  double d = a[1] + 1 / a[5], f;
  int x,y,i = 0;

  axis.el(0) = cos(a[2])*cos(a[3])*cos(a[4]) - sin(a[2])*sin(a[4]);
  axis.el(1) = sin(a[2])*cos(a[3])*cos(a[4]) + cos(a[2])*sin(a[4]);
  axis.el(2) = -sin(a[3])*cos(a[4]);

  // cout << "\nAdjust axis:" << axis;
    
  n.el(0) = cos(a[2])*sin(a[3]);
  n.el(1) = sin(a[2])*sin(a[3]);
  n.el(2) = cos(a[3]);

  n.multiply_col(0,d);

  // cout << "\nAdjust vect:" << n;
   
  for(x = reg.minx; x <= reg.maxx; x++)
    for(y = reg.miny; y <= reg.maxy; y++)
      if (reg.included(x, y))
      { p = reg.get_point(x,y);
        pt.el(0) = p.x - n.el(0) - a[11];
        pt.el(1) = p.y - n.el(1) - a[12];
        pt.el(2) = p.z - n.el(2) - a[13];
        f = axis.inner(pt);
        if (!i++) 
	{ zmax = zmin = f;
	  } else
	  { if (f < zmin) zmin = f;
              else if (f > zmax) zmax = f;
	    }
	}  
  }

void cylinder::fprint(FILE *f)
{ fprintf(f,"5 %lf %lf %lf %lf %lf %lf %lf %lf",a[1],a[2],a[3],a[4],a[5],a[11],a[12],a[13]);
  }

  
void cylinder::parameters(char **name, double *value)
{ int i;
  sprintf(name[0],"rho");
  sprintf(name[1],"phi");
  sprintf(name[2],"theta");
  sprintf(name[3],"alpha");
  sprintf(name[4],"k");
  for (i = 0; i < 5; i++)
    value[i] = a[i+1];
  value[5] = a[11];
  value[6] = a[12];
  value[7] = a[13];
  value[8] = zmin;
  value[9] = zmax;
  }

void cylinder::set_parameters(double *value)
{ int i;
  for (i = 0; i < 5; i++)
    a[i+1] = value[i];

  a[11] = value[5];
  a[12] = value[6]; 
  a[13] = value[7];
  zmin = value[8];
  zmax = value[9];
  }

struct point cylinder::project(struct point &p)
{ vector3D axis,n,v,x,t,mr;
  struct point p1;
  double d; 

  axis.el(0) = cos(a[2])*cos(a[3])*cos(a[4]) - sin(a[2])*sin(a[4]);
  axis.el(1) = sin(a[2])*cos(a[3])*cos(a[4]) + cos(a[2])*sin(a[4]);
  axis.el(2) = -sin(a[3])*cos(a[4]);
  n.el(0) = cos(a[2])*sin(a[3]);
  n.el(1) = sin(a[2])*sin(a[3]);
  n.el(2) = cos(a[3]);
  
  v.el(0) = p.x - a[11];
  v.el(1) = p.y - a[12];
  v.el(2) = p.z - a[13];

  d = a[1] + 1 / a[5];
  n.multiply_col(0,d);

  x = v - n;          // point on the axis
  t = axis;
  t.multiply_col(0,x.inner(axis));
  t += n;
  mr = v - t;
  d = fabs(1/a[5]) - mr.norm();
  mr = mr.normalize();
  mr.multiply_col(0,d);

  p1.x = p.x + mr.el(0);
  p1.y = p.y + mr.el(1);
  p1.z = p.z + mr.el(2);
  
  return(p1);
  }

void cylinder::rotate(symatrix& A)
{ class vector axis(4),n(4),o(4),b(2),c(2);
  symatrix B(2),C(2);

  axis.el(0) = cos(a[2])*cos(a[3])*cos(a[4]) - sin(a[2])*sin(a[4]);
  axis.el(1) = sin(a[2])*cos(a[3])*cos(a[4]) + cos(a[2])*sin(a[4]);
  axis.el(2) = -sin(a[3])*cos(a[4]);
  axis.el(3) = 1.0;

  n.el(0) = cos(a[2])*sin(a[3]);
  n.el(1) = sin(a[2])*sin(a[3]);
  n.el(2) = cos(a[3]);
  n.el(3) = 1.0;

  o.el(0) = a[11];
  o.el(1) = a[12];
  o.el(2) = a[13];
  o.el(3) = 1.0;

  n = A * n;
  axis = A * axis;
  o = A * o;

  a[3] = acos(n.el(2));
  a[2] = acos(n.el(0) / sin(a[3]));
  if (n.el(0) * cos(a[2]) * sin(a[3]) < 0.0) a[3] = -a[3];
    else if (n.el(1) * sin(a[2]) * sin(a[3]) < 0.0) a[2] = -a[2];

/*  B.el(0,0) = cos(a[2])*cos(a[3]);
  B.el(0,1) = -sin(a[2]);
  B.el(1,0) = sin(a[2])*cos(a[3]);
  B.el(1,1) = cos(a[2]);
  C = B.inverse(i);
  b.el(0) = axis.el(0);
  b.el(1) = axis.el(1);
  c = C * b;
  a[4] = acos(c.el(0)); */

  a[4] = acos(-axis.el(2) / sin(a[3]));
  if (fabs(axis.el(0) - (cos(a[2])*cos(a[3])*cos(a[4]) - sin(a[2])*sin(a[4]))) > 1e-6) a[4] = -a[4];
 
  

  a[11] = o.el(0);
  a[12] = o.el(1);
  a[13] = o.el(2);

  }

void cylinder::translate(vector3D &v)
{ a[11] += v.el(0);
  a[12] += v.el(1);
  a[13] += v.el(2);
  
  }


void cylinder::print()
{ vector3D axis,n;

  axis.el(0) = cos(a[2])*cos(a[3])*cos(a[4]) - sin(a[2])*sin(a[4]);
  axis.el(1) = sin(a[2])*cos(a[3])*cos(a[4]) + cos(a[2])*sin(a[4]);
  axis.el(2) = -sin(a[3])*cos(a[4]);
  n.el(0) = cos(a[2])*sin(a[3]);
  n.el(1) = sin(a[2])*sin(a[3]);
  n.el(2) = cos(a[3]);

  printf("Cylinder:  axis:(%lf,%lf,%lf)  n:(%lf,%lf,%lf)  radius:%lf\n",
               axis.el(0), axis.el(1), axis.el(2), n.el(0), n.el(1), n.el(2), 1/a[5]);
}

int cylinder::model_ok()
{ return(fabs(1/a[5]) < max_radius);
  }

void cylinder::get_normal_vector(double &x, double &y, double &z)
{ x = cos(a[2])*cos(a[3])*cos(a[4]) - sin(a[2])*sin(a[4]);
  y = sin(a[2])*cos(a[3])*cos(a[4]) + cos(a[2])*sin(a[4]);
  z = -sin(a[3])*cos(a[4]);  
  }

double cylinder::max_radius = 1e20;

