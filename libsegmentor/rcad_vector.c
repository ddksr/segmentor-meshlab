/* vector.c ---- vector manipulation routines */

#include <math.h>

double dot ( double *a, double *b)

{  double ans;

   ans = a[1]*b[1] + a[2]*b[2] + a[3]*b[3];

   return ans;
}

double modsq(double *a)

{  double ans;

   ans =  a[1]*a[1] + a[2]*a[2] + a[3]*a[3];

    return ans;
}


void normalise(double *a, double *norm)


{  double n;

    n = sqrt(a[1]*a[1] + a[2]*a[2] + a[3]*a[3]);

    norm[1] = a[1]/n;
    norm[2] = a[2]/n;
    norm[3] = a[3]/n;
}



void cross_prod( double *a, double *b, double *cross)

{   cross[1] = a[2]*b[3] - a[3]*b[2];
    cross[2] = a[3]*b[1] - a[1]*b[3];
    cross[3] = a[1]*b[2] - a[2]*b[1];
}


void vector_mult_sc(double *v, double s)

/* multiply vector by a scalar */


{   int i;
   
   for (i=1;i<4;++i)
       v[i] = v[i]*s;

}

double dist(double x1, double x2, double y1, double y2,double z1, double z2)


{  double ans,dx,dy,dz;

     dx = x1 -x2;
     dy = y1 -y2;
     dz = z1 - z2;

     ans = sqrt( dx*dx + dy*dy + dz*dz);
     
   return(ans);
}
