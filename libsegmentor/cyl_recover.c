 #include <math.h>
#include <sys/types.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

#define TINY 1.0e-20;

#define RIDICULOUS_NUMBER -1000000000



#define SQR( x ) ((x)*(x))

#define N 20000

#define PARAM_SIZE 5


void djacobi(), deigsrt(), inverse(), lubksb(), ludcmp(), nerror(),free_dvector();

double *dvector(),getval(), drand48(), getsmval(), **dmatrix();

void  o_mrqmin();

double getval(), drand48();

int *ivector();

void free_dvector(), free_ivector(), free_dmatrix();

double dot(), modsq(), dist();

void normalise(double *a, double *norm);
void cross_prod(double *a, double *norm, double *cross_norm);


 int cyldx(x,y,z,a,f,dfda,na)
 
 double x,y,z,a[],*f,dfda[];
int na;
 
 {   double *p  = dvector(1,3);
     double *acyl  = dvector(1,3);
     double *r0  = dvector(1,3);
     double *ry0bar  = dvector(1,3);
     double *ry0  = dvector(1,3);
     double *rv0  = dvector(1,3);
     double *v = dvector(1,3);
   
        double pdota, pdotr0;

     if (na) na = na; /* to avoid stupid compiler warnings */

        p[1] = x;
        p[2] = y;
        p[3] = z;

        r0[1] = cos(a[2])*sin(a[3]);
        r0[2] = sin(a[2])*sin(a[3]);
        r0[3] = cos(a[3]);

        ry0bar[1] = -sin(a[2]);
        ry0bar[2] = cos(a[2]);
        ry0bar[3] = 0;

        ry0[1] = -sin(a[2])*sin(a[3]);
        ry0[2] = cos(a[2])*sin(a[3]);
        ry0[3] = 0;

        rv0[1] = cos(a[2])*cos(a[3]);
        rv0[2] = sin(a[2])*cos(a[3]);
        rv0[3] = -sin(a[3]);


        acyl[1] = cos(a[2])*cos(a[3])*cos(a[4])  - sin(a[2])*sin(a[4]);
        acyl[2] = sin(a[2])*cos(a[3])*cos(a[4])  + cos(a[2])*sin(a[4]);
        acyl[3] = -sin(a[3])*cos(a[4]);


        pdota = dot(p,acyl);
        pdotr0 = dot(p,r0);
 
      /* calculate fit on input params - repaired by Dave 16.2.1998 */

       *f = a[5]*(modsq(p) - 2*a[1]*pdotr0 - pdota*pdota + a[1]*a[1])/2.0 + a[1] - pdotr0;
        
      /* d/drho */   
     dfda[1] =  a[5]*(a[1] -    pdotr0) + 1;

   /* d/dphi */
      v[1] = -sin(a[2])*cos(a[3])*cos(a[4])  - cos(a[2])*sin(a[4]);
      v[2] = cos(a[2])*cos(a[3])*cos(a[4])  + sin(a[2])*sin(a[4]);
      v[3] = 0;

     dfda[2] = -a[5]*(a[1]*dot(p,ry0bar) *sin(a[3]) + pdota*dot(p,v)) - dot(p,ry0);

    /* d/dtheta */

     dfda[3] = a[5]*( pdota*pdotr0*cos(a[4]) - a[1]*dot(p,rv0)) - dot(p,rv0);

    /* d/dalpha */

      v[1] = rv0[1]*sin(a[4]) - ry0bar[1]*cos(a[4]);
      v[2] = rv0[2]*sin(a[4]) - ry0bar[2]*cos(a[4]);
      v[3] = rv0[3]*sin(a[4]) - ry0bar[3]*cos(a[4]);


     dfda[4] =  a[5]*pdota*dot(p,v); 

    /* d/dk */

     dfda[5] =  0.5*(modsq(p) - 2*a[1]*pdotr0 - pdota*pdota + a[1]*a[1]);

    free_dvector(p,1,3);
    free_dvector(acyl,1,3);

    free_dvector(v,1,3);
    free_dvector(r0,1,3);
    free_dvector(ry0,1,3);
    free_dvector(ry0bar,1,3);
    free_dvector(rv0,1,3);

     return 1;
 }


int cyl_fit(double ** x, double *y, double **norm,  int ndata, double *a,  double *chisq)
{
  int i,k,index;
  double dx,d;
  int n = 0, ncyl = 0, firstflag = 0;
  
  double    xp,yp,zp;
  
 
  double    *acyl= (double *) dvector(1,3); /* average vector directions of Vector along cyl axis */
  double  *norm1 = (double *) dvector(1,3), *norm2 = (double *) dvector(1,3), *cross_norm = (double *) dvector(1,3);
   double *first_cross = (double *) dvector(1,3);
 double    *r0= (double *) dvector(1,3);
 double *q= (double *) dvector(1,3);  /* cylinder coord sys origin */

     
     
   double  *b;

    time_t t1; /* random number seed */
    
  if (chisq) chisq = chisq;   /* to avoid stupid compiler warnings */

/* numerical recipes matrices -- weird things */
/* One day unify code for all matrices equal */

         
    /* non linear stuff */
      
      
      double **list, *F, *sig, *rv0, *ry0bar;
      int *lista;
      double alamda, newchisq, oldchisq;

      
      double **covar = (double **) dmatrix(1,PARAM_SIZE,1,PARAM_SIZE),**alpha = (double **) dmatrix(1,PARAM_SIZE,1,PARAM_SIZE); 
           

    b = (double *) dvector(1,PARAM_SIZE);
  
  
  /* non linear est */
           
           list = (double **) dmatrix(1,N,1,3);
           F = (double *) dvector(1,N);
           sig = (double *) dvector(1,N);
           n=0;

                         
            rv0  =  (double *) dvector(1,3);
            ry0bar   =  (double *) dvector(1,3);

 

if ( a[1] == RIDICULOUS_NUMBER )
  { /* NEW REGION SO CONSTRUCT INITIAL ESTIMATES */

   /* set random number seed for below ops */

   (void) time(&t1);
  
   srand48((long) t1); /* use time in seconds to set seed */

   acyl[1] = acyl[2] = acyl[3] = 0.0;

  for(i=1; i<=ndata; i++) {
    xp = x[i][1];
    yp = x[i][2];
    zp = y[i];
    
    F[n+1] = 0.0;
    sig[n+1] = 1.0;
           
    norm1[1] = norm[i][1];
    norm1[2] = norm[i][2];
    norm1[3] = norm[i][3];

     index = (int ) ndata*drand48() + 1;

    norm2[1] = norm[index][1];
    norm2[2] = norm[index][2];
    norm2[3] = norm[index][3];


     cross_prod(norm1,norm2,cross_norm);

     if ( modsq(cross_norm) >  0.0001 )
        { 
            normalise(cross_norm, cross_norm);

             if (!firstflag) {
                firstflag = 1;
                first_cross[1] = cross_norm[1];
                first_cross[2] = cross_norm[2];
                first_cross[3] = cross_norm[3];
             }
             else
                  if ( dot(first_cross,cross_norm) < 0.0)
              {  cross_norm[1] *= -1.0;
                 cross_norm[2] *= -1.0;
                 cross_norm[3] *= -1.0;
            }   


     acyl[1] += cross_norm[1];
     acyl[2] += cross_norm[2];
     acyl[3] += cross_norm[3];
   
   ++ncyl;

  }
           
    ++n;
           
           if (n == N)
            { printf("Warning TOO MANY POINTS FOR STORAGE\n\n\t\tStopping point inclusion but proceeding with estimate\n");
              n= n-1;
              break;
            }
          }

    acyl[1] /= ncyl;
    acyl[2] /= ncyl;
    acyl[3] /= ncyl;


    normalise(acyl,acyl);


      /* gather estimates for least squares */



     /* make r0 normal at first point */

      r0[1] = norm[1][1];
      r0[2] = norm[1][2];
      r0[3] = norm[1][3];

     /* make rho = 1 */

      b[1] = -1.0;

      /* compute estimates of angles */ 

     /* phi = atan2(r0[2]/ro[1]); */

       b[2] = atan2(r0[2],r0[1]);


     /* theta = acos(r0[3]) */

         b[3] = acos(r0[3]);


    /*  alpha  = acos(acyl*rv0) */

        rv0[1] = cos(b[2])*cos(b[3]);
        rv0[2] = sin(b[2])*cos(b[3]);
        rv0[3] = -sin(b[3]);

        ry0bar[1] = -sin(b[2]);
        ry0bar[2] = cos(b[2]);
        ry0bar[3] = 0;


        /* normalise(rv0,rv0);  */

         dx = dot(acyl,rv0);
         d = dot(acyl,ry0bar);

         b[4] = atan2(d,dx);

   
    /* Make origin of cylinder coordinate system a point 1 "unit" off 
         1st cylinder point along normal at the point */

      q[1] = x[1][1] -  b[1]* r0[1];
      q[2] = x[1][2] -  b[1]* r0[2];
      q[3] = y[1] -  b[1]* r0[3];


    /* finally compute k */

    xp = yp = 0.0;  /* use xp,yp to gather summs for least sq fit of k (b[5]) */

   for(i=1; i<=ndata; i++) {
       /* use norm1 to store point  */
       norm1[1] = x[i][1] ;
       norm1[2] = x[i][2] ;
       norm1[3] = y[i];

       /* use norm2 to store vector pi - r0 - q */

      norm2[1] = norm1[1] - b[1]*r0[1] - q[1];
      norm2[2] = norm1[2] - b[1]*r0[2] - q[2];
      norm2[3] = norm1[3] - b[1]*r0[3] - q[3];

    /* use cross_norm to store norm2 cross acyl */

   cross_prod(norm2,acyl,cross_norm);

   dx = modsq(cross_norm);

   xp  += dot(norm2,r0)*dx;

   yp += dx*dx;
 

   /* translate points to new coord system */

    list[i][1] = norm1[1] - q[1];
    list[i][2] = norm1[2] - q[2];
    list[i][3]=  norm1[3] - q[3];

  }

  b[5] = 2*xp/yp;  

}

else

{  /* copy current fit parameters to b and q */

         b[1] = a[1];
         b[2] = a[2];
         b[3] = a[3];
         b[4] = a[4];
         b[5] = a[5];

         q[1] = a[11];
         q[2] = a[12];
         q[3] = a[13];


         for(i=1; i<=ndata; i++) {
                     /* set NR variables */
                    F[n+1] = 0.0;
                    sig[n+1] = 1.0;

             
              
             /* translate points to now coord system */

             list[i][1] = x[i][1] - q[1];
             list[i][2] = x[i][2] - q[2];
             list[i][3]=  y[i] - q[3];

             ++n;

        }
}





    
                    
         
           lista = (int *) ivector(1,PARAM_SIZE);
          
           lista[1] = 1;
           lista[2] = 2;
           lista[3] = 3;
           lista[4] = 4;
           lista[5] = 5;
           
           alamda = -1.0;
           newchisq = oldchisq = 1.0;

           /* printf("\nInitial estimates:rho=%8.4lf phi=%8.4lf theta=%8.4lf alpha=%8.4lf k=%8.4lf\n",
             b[1],b[2],b[3],b[4],b[5]); */
           
           o_mrqmin(list, F, sig, n,b,PARAM_SIZE, lista,PARAM_SIZE, covar, alpha, &newchisq,&alamda,cyldx); 
           
            k= 0;
           
           while (  k++ < 25 )
            { 
              oldchisq = newchisq;
              o_mrqmin(list, F, sig, n,b,PARAM_SIZE, lista,PARAM_SIZE, covar, alpha, &newchisq,&alamda,cyldx);
              
                           
            }
            
           alamda = 0.0;
           o_mrqmin(list, F, sig, n,b,PARAM_SIZE, lista,PARAM_SIZE, covar, alpha, &newchisq,&alamda,cyldx);
           
                     
                 
           a[1] = b[1];
           a[2] = b[2];
           a[3] = b[3]; 
           a[4] = b[4]; 
           a[5] = b[5];

           a[11] = q[1];
           a[12] = q[2];
           a[13] = q[3];

           
         

      

 
    free_dvector(acyl,1,3);
    free_dvector(norm1,1,3);
    free_dvector(norm2,1,3);
    free_dvector(cross_norm,1,3);

    free_dvector(r0,1,3);
    free_dvector(q,1,3);

    free_dvector(first_cross,1,3);



  
  free_dmatrix(list,1,N,1,3);
  free_dmatrix(covar,1,PARAM_SIZE,1,PARAM_SIZE);
  free_dmatrix(alpha,1,PARAM_SIZE,1,PARAM_SIZE);
  
  free_ivector(lista,1,PARAM_SIZE);
  free_dvector(b,1,PARAM_SIZE);
  free_dvector(F,1,N);
  free_dvector(sig,1,N);

    free_dvector(rv0,1,3);
    free_dvector(ry0bar,1,3);

   





return (2);
}



/* NEEDS TO COPE WITH A LIST NOT A SQUARE ARRAY

double getval(y,i,j,offset)
 
   int i,j;
   double offset;
   
   {     
     if ( buffer[i+j*(int)width] == 0 )
       return(0.0);
    else
      return (double) ((double)buffer[i+j*(int)width] - offset);
   }

*/




   


    





















