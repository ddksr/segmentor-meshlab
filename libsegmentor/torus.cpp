// torus.C

#include "torus.h"
#include "sline.h"
#include "pplane.h"

#include <stdlib.h>
#include <iostream>
#include <math.h>


extern "C" {

  #include "rcad_vector.h"
  #include "rcad_ndm_ut.h"
  #include "rcad_circurv.h"

  #include "num_anal.h"

void o_mrqmin(double **x, double *y,double *sig,int ndata,double *a,int ma,int *lista,int mfit,double **covar,double **alpha,double *chisq,double *alamda, int (*funcs)(double,double ,double, double *,double *,double *,int));

int cyl_fit(double **list_x, double *list_y, double **norm, int no, double *param, double *chisq);
}

#define TRUE 1
#define FALSE 0

#define PARAM_SIZE 8

#define KTHRESH  0.1  // tolerance of variance of k curvature estimate

#define CYL_THRESH 0.0001  // threshold for s curvature estimate 
                                              // below which degenerate to cylinder fitting
#define TCTRESH 10000.0

#define RIDICULOUS_NUMBER -100000000000.0

// double torus::m_error = max_description_error;   // provide maximum point distance and description error
// double torus::m_dist = max_point_distance;    // for all objects of the class

extern double gmax_description_error, gmax_point_distance;


torus::torus(region &reg) {
 
int Max = reg.point_count()+5; 


  // interface to C code and Numerical Recipes

  double **list_x = dmatrix(1,Max,1,2);
  double *list_y  = dvector(1,Max);
  double *param = dvector(1,13);
  double chisq;
 
  double **list_norm = dmatrix(1,Max,1,3);


  int x, y, min_x,min_y,max_x,max_y, no;
  struct point p;

  no = 0;
  

  min_x = reg.minx;
  min_y = reg.miny;
  max_x = reg.maxx;
  max_y = reg.maxy;
  
  
       
 
   for(x = min_x; x <= max_x; x++)
   for(y = min_y; y <= max_y; y++)
     if (reg.included(x, y))
      if (no < Max) { 
	  no++;
 	  p = reg.get_point(x, y);
 	  list_x[no][1] = p.x;
 	  list_x[no][2] = p.y;
	  list_y[no] = p.z; 

          list_norm[no][1] = p.nx;
 	  list_norm[no][2] = p.ny;
	  list_norm[no][3] = p.nz;

	  
	  
	 
 	}	
 	else
  	  std::cout << "Over the limit of Torus points";

  param[1] = RIDICULOUS_NUMBER - 1.0;

int xspan = max_x - min_x;
int yspan = max_y - min_y;

if( torus_fit(list_x, list_y, list_norm, no, xspan, yspan, param, &chisq, reg) )
   { /* Gabor rep */
     a[1] = param[1];
     a[2] = param[2];
     a[3] = param[3];
     a[4] = param[4];
     a[5] = param[5];
     a[6] = param[6];
     a[7] = param[7];
     a[8] = param[8];



     a[11] = param[11];
     a[12] = param[12];
     a[13] = param[13];
}

  
  
  free_dmatrix(list_norm,1,Max,1,3);
  free_dvector(param,1,13);
  free_dvector(list_y,1,Max);
  free_dmatrix(list_x,1,Max,1,2);
}

torus::torus(region &reg, torus &init) {
 
int Max = reg.point_count()+1;  


  // interface to C code and Numerical Recipes

  double **list_x = dmatrix(1,Max,1,2);
  double *list_y  = dvector(1,Max);
  double *param = dvector(1,13);
  double chisq;

   double **list_norm  = dmatrix(1,Max,1,3);


  int x, y, min_x,min_y,max_x,max_y, no;
  struct point p;

  no = 0;
  
  min_x = reg.minx;
  min_y = reg.miny;
  max_x = reg.maxx;
  max_y = reg.maxy;
  

     param[1] = init.a[1];
     param[2] = init.a[2];
     param[3] = init.a[3];
     param[4] = init.a[4];
     param[5] = init.a[5];
     param[6] = init.a[6];
     param[7] = init.a[7];
     param[8] = init.a[8];




     param[11] = init.a[11];
     param[12] = init.a[12];
     param[13] = init.a[13];


   

     
 
  for(x = min_x; x <= max_x; x++)
   for(y = min_y; y <= max_y; y++)
      if (reg.included(x, y))
	if (no < Max) { 
	  no++;
	  p = reg.get_point(x, y);
	  list_x[no][1] = p.x;
	  list_x[no][2] = p.y;
	  list_y[no] = p.z;

          list_norm[no][1] = p.nx;
 	  list_norm[no][2] = p.ny;
	  list_norm[no][3] = p.nz;

	 
	}	
	else
	  std::cout << "Over the limit of Cylinder points";


          
int xspan = max_x - min_x;
int yspan = max_y - min_y;

if( torus_fit(list_x, list_y, list_norm, no, xspan, yspan, param, &chisq, reg) )
 
   { /* Gabor rep */
     a[1] = param[1];
     a[2] = param[2];
     a[3] = param[3];
     a[4] = param[4];
     a[5] = param[5];
     a[6] = param[6];
     a[7] = param[7];
     a[8] = param[8];



     a[11] = param[11];
     a[12] = param[12];
     a[13] = param[13];

}

     
 
  free_dmatrix(list_x,1,Max,1,2);
  free_dvector(list_y,1,Max);
  free_dvector(param,1,13);
  free_dmatrix(list_norm,1,Max,1,3);
}

// Euclidian distance - implemented by Bojan Kverh 9.6.1998

double torus::abs_signed_distance(struct point &p)
{ double sinphi = sin(a[2]), cosphi = cos(a[2]);
  double sintheta = sin(a[3]), costheta = cos(a[3]);
  double sinsigma = sin(a[4]), cossigma = cos(a[4]);
  double sintau = sin(a[5]), costau = cos(a[5]);
  double d,x1,x2;

  vector3D n,pt,axis,n1,n2,p1,p2,ncrossa,p2crossa;

  n.el(0) = cosphi*sintheta;
  n.el(1) = sinphi*sintheta;
  n.el(2) = costheta;

  pt.el(0) = p.x - a[11];
  pt.el(1) = p.y - a[12];
  pt.el(2) = p.z - a[13];

  axis.el(0) = cossigma*sintau;
  axis.el(1) = sinsigma*sintau;
  axis.el(2) = costau;

  n2 = n1 = n;
  n1.multiply_col(0,a[1]+1/a[6]);
  n2.multiply_col(0,a[1]+1/a[7]);
   
  p1 = pt - n1;
  p2 = pt - n2;

  ncrossa = n.outer(axis);
  p2crossa = p2.outer(axis);

  x1 = p1.inner(axis);
  x2 = p2crossa.norm()-a[8]*fabs(1/a[7]-1/a[6])*ncrossa.norm();

  // std::cout << "\nDist:" << x1 << " " << x2 << " ";

  d = sqrt(x1*x1+x2*x2)-1/fabs(a[6]);  

  return(d);
  }

/*
double torus::distance(struct point &v)
{  double r_est;
   double modsqncrossa;  

    double sinphi = sin(a[2]), cosphi = cos(a[2]);
    double sintheta = sin(a[3]), costheta = cos(a[3]);
    double sinsigma = sin(a[4]), cossigma = cos(a[4]);
    double sintau = sin(a[5]), costau = cos(a[5]);

      
   double *p  = dvector(1,3);
   double *phat  = dvector(1,3);
   double *atorus  = dvector(1,3);
   double *n  = dvector(1,3);
   double *ncrossa = dvector(1,3);
   double *phatminusndivs = dvector(1,3);
   double *phatminusndivscrossa = dvector(1,3);



   double phatdota, phatdotn, pdota, pdotn, pdotncrossa, ndota,phatminusndivscrossadotncrossa;

    double sphere_dist,delta_dist;
    
     int sign;


    p[1] = v.x - a[11];
    p[2] = v.y - a[12];
    p[3] = v.z - a[13];

    // calculate n

    n[1] = cosphi*sintheta;
    n[2] = sinphi*sintheta;
    n[3] = costheta;

    // phat - p - rho*n

     phat[1] = p[1] - a[1]*n[1];
     phat[2] = p[2] - a[1]*n[2];
     phat[3] = p[3] - a[1]*n[3];

     // compute phatminusndivs

    phatminusndivs[1] = phat[1] - n[1]/a[7];
    phatminusndivs[2] = phat[2] - n[2]/a[7];
    phatminusndivs[3] = phat[3] - n[3]/a[7];


    // compute atorus

     atorus[1] = cossigma*sintau;
     atorus[2] = sinsigma*sintau;
     atorus[3] = costau;

    
     // compute cross  products

    cross_prod(n,atorus,ncrossa);
    modsqncrossa = modsq(ncrossa);

    cross_prod(phatminusndivs,atorus, phatminusndivscrossa); 


      // compute dot products

     phatdota = dot(phat,atorus);
     phatdotn = dot(phat, n);
     ndota = dot(n,atorus);
     pdotn = dot(p, n);
     pdota = dot(p, atorus);
     pdotncrossa = dot(p, ncrossa);


     phatminusndivscrossadotncrossa = dot(phatminusndivscrossa,ncrossa);


    sphere_dist = 0.5*a[6]*modsq(phat) - phatdotn;

    // work out sign of k*k/s -k: use r_est to store value of

   r_est = a[6]*a[6]/a[7] - a[6];
if ( fabs(r_est) < CYL_THRESH)   // use thresh just to check for near zero here
     sign = 0;
else
  if ( r_est < 0)
    sign = -1;
  else
    sign = 1;

    
    
 if  ( fabs(a[7]) > CYL_THRESH) 

   delta_dist = (a[6]/a[7] - 1.0)*( a[8]*sign* sqrt(modsq(phatminusndivscrossa))*sqrt(modsqncrossa) + phatminusndivscrossadotncrossa); 

    
    else // degenerate to cylinder distance to avoid div by 0;
     {   // using cylinder formulation use atorus for cylinder axis

        pdota = dot(p, atorus);

        atorus[1] = cos(a[2])*cos(a[3])*cos(a[4])  - sin(a[2])*sin(a[4]);
        atorus[2] = sin(a[2])*cos(a[3])*cos(a[4])  + cos(a[2])*sin(a[4]);
        atorus[3] = -sin(a[3])*cos(a[4]);


       delta_dist = a[6]*(modsq(p) - 2*a[1]*pdotn - pdota*pdota + a[1]*a[1])/2.0 + a[1] - pdotn;

     }

    r_est = sphere_dist - delta_dist;

    // tidy up memory of nr utils

    free_dvector(p,1,3);
    free_dvector(phat,1,3);
    free_dvector(n,1,3);
    free_dvector(atorus,1,3);
    free_dvector(ncrossa,1,3);
    free_dvector(phatminusndivs,1,3);
    free_dvector(phatminusndivscrossa,1,3);

    return fabs(r_est);
       
 }
  */ 


void torus::print()
{ 

if ( a[8] != 0) // regular torus
{
 
  std::cout << "Torus:\n\t\t  rho =  " << a[1] << ", phi = " << a[2] << ",  theta = " << a[3] << "\n\t\t sigma =  " << a[4] << ", tau = " << a[5] << ",  k= " << a[6] << ",  s= " << a[7] << ",  epsilon= " << a[8]<< '\n' << ",  Origin= " << a[11] << ", " << a[12] << ", " << a[13] << '\n';
 }
else
{  std::cout << "Torus(currently denerated to a Cylinder):\n\t\t  rho =  " << a[1] << ", phi = " << a[2] << ",  theta = " << a[3] << "\n\t\t alpha =  " << a[4] << ",  k= " << a[6] << '\n' << ",  Origin= " << a[11] << ", " <<
 a[12] << ", " << a[13] << '\n';
}

}





void torus::save(FILE *f) { 

  fwrite(a,14,sizeof(double),f);
  
}

torus::torus(FILE *f) {
 
    
  fread(a,14,sizeof(double),f);
 

}

model* torus::improve(region& orig, region& added)
{ region r = orig | added;
  return(new torus(r));
  }

model* torus::improve(model *m, region& orig, region& added)
{ region r = orig | added;
  torus *t;
  t = new torus(r,*(torus *)m);

  if (!t->torus_ok())
  {
    *t = *(torus *)m;    
    }

  return(t);
  
  }

int torus::torus_ok()
{ return(fabs(a[6]) < TCTRESH && fabs(a[7]) < TCTRESH);
  }


void sphdx(double param[],  double modsqp, double pdotn, double pdotdndphi, double pdotdndtheta, double dfda[])


{
// d/drho

dfda[1] = param[6]*(param[1] - pdotn) +1;

// d/dphi

dfda[2] = pdotdndphi*(-1.0*param[6]*param[1] - 1.0);

// d/dtheta

dfda[3] = pdotdndtheta*(-1.0*param[6]*param[1] - 1.0);


// d/dk

dfda[6] = 0.5*(modsqp - 2.0*param[1]*pdotn + param[1]*param[1]);

}


void deltadx(double param[], double modsqncrossa, double phatminusndivscrossa[],  double ndota, double ncrossa[], double dndphicrossa[], double dndthetacrossa[], double phatminusndivscrossdadsigma[], double phatminusndivscrossdadtau[], double ncrossdadsigma[], double ncrossdadtau[], double dfda[])


{   double *v=dvector(1,3);
    double *u=dvector(1,3);

    double vdotu, ksminus1, rhoplus1s;

    double modphatminusndivscrossa = sqrt(modsq(phatminusndivscrossa));

    double modncrossa = sqrt(modsqncrossa);

    int sign,signk;

    double epssign; // epsilon * sign(k*k/s  -k)

     double dndphicrossadotv, dndphicrossadotu, dndthetacrossadotv,dndthetacrossadotu;

     double  phatminusndivscrossadotv, phatminusndivscrossadotu;

     double phatminusndivscrossdadsigmadotu, phatminusndivscrossdadtaudotu;
       
      double ncrossdadsigmadotv, ncrossdadsigmadotu;

       double ncrossdadtaudotu, ncrossdadtaudotv;

    
    v[1] = phatminusndivscrossa[1]/modphatminusndivscrossa;
    v[2] = phatminusndivscrossa[2]/modphatminusndivscrossa;
    v[3] = phatminusndivscrossa[3]/modphatminusndivscrossa;


    u[1] = ncrossa[1]/modncrossa;
    u[2] = ncrossa[2]/modncrossa;
    u[3] = ncrossa[3]/modncrossa;


    ksminus1 = param[6]/param[7] - 1.0;

    rhoplus1s = param[1] + 1.0/param[7];

    vdotu = dot(v,u);

    dndphicrossadotv = dot(dndphicrossa,v);
    dndphicrossadotu = dot(dndphicrossa,u);
    dndthetacrossadotv = dot(dndthetacrossa,v);
    dndthetacrossadotu = dot(dndthetacrossa,u);

    phatminusndivscrossadotv = dot(phatminusndivscrossa,v);
    phatminusndivscrossadotu = dot(phatminusndivscrossa,u);

    phatminusndivscrossdadsigmadotu = dot(phatminusndivscrossdadsigma,u);
    phatminusndivscrossdadtaudotu = dot(phatminusndivscrossdadtau,u);


    ncrossdadsigmadotv =  dot(ncrossdadsigma,v);

    ncrossdadsigmadotu =  dot(ncrossdadsigma,u);

    ncrossdadtaudotu =  dot(ncrossdadtau,u);

    ncrossdadtaudotv =  dot(ncrossdadtau,v);




     

 // work out sign of k*k/s -k: use 
    
  if ( (param[6]*param[6]/param[7] - param[6]) < 0)
    sign = -1;
  else
    sign = 1;

epssign = sign*param[8];


 // signk

if ( (param[6]) < 0)
    signk = -1;
  else
    signk = 1;

   

// d/drho

   dfda[1] = -1.0*modsqncrossa*ksminus1*(epssign*vdotu + 1);

// d/dphi

    dfda[2] = -1.0*rhoplus1s*ksminus1*(epssign*modphatminusndivscrossa + modncrossa)*(epssign*dndphicrossadotv + dndphicrossadotu);

// d/dtheta

       dfda[3] = -1.0*rhoplus1s*ksminus1*(epssign*modphatminusndivscrossa + modncrossa)*(epssign*dndthetacrossadotv + dndthetacrossadotu);

// d/sigma

    dfda[4] = ksminus1*(epssign*(phatminusndivscrossdadsigmadotu*ndota + ncrossdadsigmadotu*modphatminusndivscrossa) + phatminusndivscrossdadsigmadotu*modncrossa + ncrossdadsigmadotv*modphatminusndivscrossa);

// d/dtau

    dfda[5] = ksminus1*(epssign*(phatminusndivscrossdadtaudotu*ndota + ncrossdadtaudotu*modphatminusndivscrossa) + phatminusndivscrossdadtaudotu*modncrossa + ncrossdadtaudotv*modphatminusndivscrossa);

// d/dk

    dfda[6] = modphatminusndivscrossa*modncrossa*(epssign + vdotu)/param[7];

// d/ds

    dfda[7] = modncrossa*(epssign + vdotu)*(param[8]*signk*fabs(ksminus1)*modncrossa - param[6]*modphatminusndivscrossa)/(param[7]*param[7]);



  // clean up nr memory allocs

    free_dvector(v,1,3);
    free_dvector(u,1,3);


}


  



int torusdx(double x, double y, double z, double param[], double *f, double dfda[], int na)

{  

    double *spheredfda  = dvector(1,PARAM_SIZE);
    double *deltadfda  = dvector(1,PARAM_SIZE);

   
   
    double sinphi = sin(param[2]), cosphi = cos(param[2]);
    double sintheta = sin(param[3]), costheta = cos(param[3]);
    double sinsigma = sin(param[4]), cossigma = cos(param[4]);
    double sintau = sin(param[5]), costau = cos(param[5]);

   double modsqp,modsqphat,modsqncrossa;  
      
   double *p  = dvector(1,3);
   double *p2  = dvector(1,3);
   double *phat  = dvector(1,3);
   double *atorus  = dvector(1,3);
   double *acyl = dvector(1,3);   // needed for degenerate case when s -> 0

   double *dadsigma  = dvector(1,3);
   double *dadtau  = dvector(1,3);
   double *n  = dvector(1,3);
   double *dndphi = dvector(1,3);
   double *dndtheta  = dvector(1,3);
   double *ncrossa = dvector(1,3);
   double *dndphicrossa = dvector(1,3);

   double *dndthetacrossa = dvector(1,3);

   double *ncrossdadsigma = dvector(1,3);
   double *ncrossdadtau= dvector(1,3);


   double *phatminusndivs = dvector(1,3);
   double *phatminusndivscrossa = dvector(1,3);
   double *phatminusndivscrossdadsigma = dvector(1,3);
   double *phatminusndivscrossdadtau = dvector(1,3);
 



   double pdota, pdotn,
              pdotdndphi, pdotdndtheta,
              phatdota, phatdotn, 
              phatdotdadsigma, phatdotdadtau,
              phatdotdndphi, phatdotdndtheta,
              ndota,  ndotdadsigma, ndotdadtau,
              dndphidotp,dndphidota,
              dndthetadotp,dndthetadota,
              phatminusndivscrossadotncrossa;
   
   double p2dota;

   double sphere_dist, delta_dist;

   int sign;


    if (na) na = na;  // to avoid stupid compiler warnings

    p[1] = x;
    p[2] = y;
    p[3] = z ;

    modsqp = modsq(p);

    // calculate n, dndphi, dndtheta

    n[1] = cosphi*sintheta;
    n[2] = sinphi*sintheta;
    n[3] = costheta;

    dndphi[1] = -sinphi*sintheta;
    dndphi[2] = cosphi*sintheta;
    dndphi[3] = 0;

    dndtheta[1] = cosphi*costheta;
    dndtheta[2] = sinphi*costheta;
    dndtheta[3] = -sintheta;


    // phat = p - rho*n

     phat[1] = p[1] - param[1]*n[1];
     phat[2] = p[2] - param[1]*n[2];
     phat[3] = p[3] - param[1]*n[3];

      modsqphat = modsq(phat);

     // compute phatminusndivs

    phatminusndivs[1] = phat[1] - n[1]/param[7];
    phatminusndivs[2] = phat[2] - n[2]/param[7];
    phatminusndivs[3] = phat[3] - n[3]/param[7];


    // p2 = p - 2*rho*n

     p2[1] = p[1] - 2.0*param[1]*n[1];
     p2[2] = p[2] - 2.0*param[1]*n[2];
     p2[3] = p[3] - 2.0*param[1]*n[3];

    // compute atorus, dadsigma, dadtau

     atorus[1] = cossigma*sintau;
     atorus[2] = sinsigma*sintau;
     atorus[3] = costau;

     dadsigma[1] = -sinsigma*sintau;
     dadsigma[2] = cossigma*sintau;
     dadsigma[3] = 0;

     dadtau[1] = cossigma*costau;
     dadtau[2] = sinsigma*costau;
     dadtau[3] = -sintau;


     // compute dot products
     pdota = dot(p,atorus);
     pdotn = dot(p, n);

     pdotdndphi = dot(p, dndphi);
     pdotdndtheta = dot(p, dndtheta);


     phatdota = dot(phat,atorus);
     phatdotn = dot(phat, n);
     phatdotdadsigma = dot(phat, dadsigma);
     phatdotdadtau = dot(phat, dadtau);
     phatdotdndphi = dot(phat, dndphi);
     phatdotdndtheta = dot(phat, dndtheta);

     p2dota = dot(p2,atorus);

     ndota = dot(n,atorus);
     ndotdadsigma =  dot(n,dadsigma);
     ndotdadtau =  dot(n,dadtau);

     dndphidotp = dot(dndphi,p);
     dndphidota = dot(dndphi,atorus);

     dndthetadotp = dot(dndtheta,p);
     dndthetadota = dot(dndtheta,atorus);


     // compute n cross a product

    cross_prod(n,atorus,ncrossa);
    modsqncrossa = modsq(ncrossa);


    cross_prod(phatminusndivs,atorus, phatminusndivscrossa); 

    cross_prod(dndphi,atorus,dndphicrossa);
    cross_prod(dndtheta,atorus,dndthetacrossa);


    cross_prod(phatminusndivs,dadsigma, phatminusndivscrossdadsigma); 

    cross_prod(phatminusndivs,dadtau, phatminusndivscrossdadtau); 

    cross_prod(n,dadsigma, ncrossdadsigma);
    cross_prod(n,dadtau, ncrossdadtau);


     phatminusndivscrossadotncrossa = dot(phatminusndivscrossa,ncrossa);



  

    if  ( fabs(param[7]) > CYL_THRESH ) 
        {  // Normal torus fit procedure

           // compute sphere partial derivatives

   sphdx(param, modsqp,    pdotn,  pdotdndphi,  pdotdndtheta,  spheredfda);

  //compute delta partial derivatives


   deltadx( param, modsqncrossa,  phatminusndivscrossa,   ndota, ncrossa,  dndphicrossa,  dndthetacrossa,  phatminusndivscrossdadsigma,  phatminusndivscrossdadtau,  ncrossdadsigma,  ncrossdadtau,  deltadfda);

        sphere_dist = 0.5*param[6]*modsq(phat) - phatdotn;

    // work out sign of k*k/s -k: use r_est to store value of

   *f = param[6]*param[6]/param[7] - param[6];

if ( fabs(*f) < CYL_THRESH)  // use thresh just to check for near zero
     sign = 0;
else
  if ( *f < 0)
    sign = -1;
  else
    sign = 1;

   delta_dist = (param[6]/param[7] - 1.0)*( param[8]*sign* sqrt(modsq(phatminusndivscrossa))*sqrt(modsqncrossa) + phatminusndivscrossadotncrossa); 

   
    *f = sphere_dist - delta_dist;


// compute  torus dfda 

   // d/drho

   dfda[1] = spheredfda[1] - deltadfda[1];

// d/dphi

    dfda[2] = spheredfda[2] - deltadfda[2];

// d/dtheta

       dfda[3] = spheredfda[3] - deltadfda[3];
 
// d/sigma

    dfda[4] =  -1.0*deltadfda[4];

// d/dtau

    dfda[5] = -1.0*deltadfda[5];

// d/dk

    dfda[6] = spheredfda[6] - deltadfda[6];

// d/ds

    dfda[7] = -1.0*deltadfda[7];

  }
else 

   {   // degenerate to cylinder distance to avoid div by 0;

        printf("Should  NOT reach hear: Curvature s below threshold in torusdx\n");
        return(0);

    } 
 
   // tidy up memory of nr utils

    free_dvector(spheredfda,1,PARAM_SIZE);
    free_dvector(deltadfda,1,PARAM_SIZE);

    free_dvector(p,1,3);
    free_dvector(p2,1,3);
    free_dvector(phat,1,3);
    free_dvector(atorus,1,3);
    free_dvector(acyl,1,3);

    free_dvector(dadsigma,1,3);
    free_dvector(dadtau,1,3);
    free_dvector(n,1,3);
    free_dvector(dndphi,1,3);
    free_dvector(dndtheta,1,3);
    free_dvector(ncrossa,1,3);
    free_dvector(dndphicrossa,1,3);

    free_dvector(dndthetacrossa,1,3);
    free_dvector(ncrossdadsigma,1,3);

    free_dvector(ncrossdadtau,1,3);


    free_dvector(phatminusndivs,1,3);
    free_dvector(phatminusndivscrossa,1,3);
    free_dvector(phatminusndivscrossdadsigma,1,3);
    free_dvector(phatminusndivscrossdadtau,1,3);



    return 1;
 
}

int torus::torus_fit(double **list_x, double *list_y, double **norm, int no, int xspan, int yspan,double *param, double *chisq, region &reg)

{
  int ni,nj,ki,kj,ix,iy;
  double pt[3], pt_norm[3];
  double  atorus[3];

  ix = iy = 0;
      
  region *rnew;

  int ind[9],converted[9];
  sline input_norms[5];
  sline output_axis;
  sline *r_axis=&output_axis;
  int   axnum=1;

  int i, index, j, k, index2, index_base;
  double angle, max_angle = -10000000.0;

   double *k1,*k2, temp;
    int torus_vote[2]; // index 0 collect votes for lemon torus, index 1 apple torus
    int step, step_no, x_off, y_off;
     int error;
    double sumk1, sumk1sq, vark1, sumk2, sumk2sq, vark2;   /// stats to test torus



   double *q= (double *) dvector(1,3);  /* torus coord sys origin */

 // Numerical recipe vars
     
     double **list, *F, *sig;
     double  *b = dvector(1,PARAM_SIZE);
     double  *n = dvector(1,3);
      int *lista;
      double alamda, newchisq, oldchisq;

      
      double  **covar = (double **) dmatrix(1,PARAM_SIZE,1,PARAM_SIZE),
                  **alpha = (double **) dmatrix(1,PARAM_SIZE,1,PARAM_SIZE); 


      rnew = new region(reg.theImage);

      list = (double **) dmatrix(1,no,1,3);
      F = (double *) dvector(1,no);
      sig = (double *) dvector(1,no);

       k1 = (double *) dvector(1,no);
       k2 = (double *) dvector(1,no);

     pplane projection_plane(reg,ind),*pp1;

     if (yspan) yspan = yspan;   // to avoid stupid compiler warnings.

if ( ( param[1] <= RIDICULOUS_NUMBER)  ||  (b[8] == 0) )
  { /* NEW REGION SO CONSTRUCT INITIAL ESTIMATES */
        
    // b[8] == 0 if current fit is a cylinder (degenerate case) if so try new patch and reestimate
    //  initial parameters.

  // Guess AXIS of Torus 

  //  Take 5 pts at 4 corners and center of region and compute sline of normals

    // top left normal

    index = ind[0];

    pt[0] =  list_x[index][1];
    pt[1] = list_x[index][2];
    pt[2] = list_y[index];

    for (i=0;i<3;++i)
      pt_norm[i] = norm[index][i+1];

   input_norms[0] = sline(pt, pt_norm,0); 

    // top right normal


   index = ind[2];

    pt[0] =  list_x[index][1];
    pt[1] = list_x[index][2];
    pt[2] = list_y[index];

    for (i=0;i<3;++i)
      pt_norm[i] = norm[index][i+1];


   input_norms[1] = sline(pt, pt_norm,0);
 
   // middle normal

   index =  ind[4];

    pt[0] =  list_x[index][1];
    pt[1] = list_x[index][2];
    pt[2] = list_y[index];

    for (i=0;i<3;++i)
      pt_norm[i] = norm[index][i+1];

   input_norms[2] = sline(pt, pt_norm,0); 

    // bottom left normal

   index = ind[6];

    pt[0] =  list_x[index][1];
    pt[1] = list_x[index][2];
    pt[2] = list_y[index];

    for (i=0;i<3;++i)
      pt_norm[i] = norm[index][i+1];

   input_norms[3] = sline(pt, pt_norm,0); 

   // bottom right normal

    index = ind[8];

    pt[0] =  list_x[index][1];
    pt[1] = list_x[index][2];
    pt[2] = list_y[index];

    for (i=0;i<3;++i)
      pt_norm[i] = norm[index][i+1];

   input_norms[4] = sline(pt, pt_norm,0); 

   index2 = index = guess_axes(input_norms,5,
             &r_axis, &axnum,
             1.5,
             NULL,
             10,
             1000.0,
             0.001,
             0.0005,
             1,
             NULL,
             0.5,
             NULL);

    output_axis.direction(atorus);

 // index = 1;
 // atorus[0] = atorus[2] = 0.0;
//  atorus[1] = 1.0;

 if ( (index == 0) || ( fabs(modsq(atorus-1) - 0.999999) > 0.001 ))
       {  // we have no axis
         std::cout << "ERROR: Cannot Guess Axis for this Torus - guessing in the primitive way\n";
         // return(0);

       }

 // check for axis - if no axis can be guessed return FALSE

  if (atorus[0] == 0.0 && atorus[1] == 0.0 && atorus[2] == 0.0)
    {          
      // Added by Bojan Kverh 1.6.1998    
      // Estimation of torus axis
      // If this can't be guessed by guess_axes, try this.
      
      vector3D v1,v2,v3,ax;
      double dnorm;
      int npts = 0;      

      index = 1 + rand() % no;

      ax.el(0) = 0.0;
      ax.el(1) = 0.0;
      ax.el(2) = 0.0; 

      for (i = 1; i <= no; i++)
      { index = 1 + rand() % no;  // random number
        v1.el(0) = norm[index][1]; 
        v1.el(1) = norm[index][2];
        v1.el(2) = norm[index][3];
        v2.el(0) = norm[i][1];
        v2.el(1) = norm[i][2];
        v2.el(2) = norm[i][3];
        v3 = v1.outer(v2);
        dnorm = v3.norm();
        if (dnorm > 0.1)
	{ v3.normalize();
          npts++;
          // if (ax.inner(v3) < 0.0) v3.multiply_col(0,-1);
          ax += v3;
	  }
	}
      
      atorus[0] = ax.el(0) / npts;
      atorus[1] = ax.el(1) / npts;
      atorus[2] = ax.el(2) / npts;

      // End Added by Bojan Kverh 1.6.1998 
      }


  // printf("%d:%f %f %f\n",index,atorus[0],atorus[1],atorus[2]);

  // compute  initial estimates

   // find base point choose point with largest angle of normal and axis.

       step = (int)   ( sqrt((double) no) / 2);

     index = index_base = index2;

     for (i = 1; i< no; i += step)
        {  angle =  acos( dot(norm[i],(atorus-1)) );
           if ( max_angle < angle )
              { index_base = i;  // remember point
                max_angle = angle;
              }
         }

    index = index_base;

    n[1] = norm[index][1];
    n[2] = norm[index][2];
    n[3] = norm[index][3];

 // make rho = 1

   b[1] = -1.0;

     /* Make origin of torus coordinate system a point 1 "unit" off 
          torus  base point point along normal at the point */

      q[1] = list_x[index][1] -  b[1]* n[1];
      q[2] = list_x[index][2] -  b[1]* n[2];
      q[3] = list_y[index] -  b[1]* n[3];

       

       


   
     /* compute estimates of angles */ 

     /* phi = atan2(n[2]/n[1]); */

       b[2] = atan2(n[2],n[1]);


     /* theta = acos(n[3]) */ 

         b[3] = acos(n[3]);


     /* sigma = atan2(atorus[2]/atorus[1]);  but in NORAML C ARREAY INDEX!! */

       b[4] = atan2(atorus[1],atorus[0]);


     /* tau = acos(n[3]) but in NORMAL C ARREAY INDEX!! */

       b[5] = acos(atorus[2]);

     //  To Estimate k  + s --- use ralph`s curvature estimation.
    //    SInce Torus  one curvature k should  always be negative and constant accross the patch or close to it.
     //   Scan through pacth and vote on this also compute variance in k and s vals to measure goodness. 
     
      rcad_SurfaceData pnts[9];
      int              inds[9];
    
       double normal[3];
     double dir1[3];
     double dir2[3]; 


      //  scan through patch  taking start, middle and end of each row/col as a sample curvature estimator

      // get regin centre pixel and compute span

     
      x_off = xspan/2;
      y_off = xspan/2;
      
      for (j=0;j<9;++j)
         inds[j] = j;

     for  (i = 1; i <= no; i += step)
        { 
          /* // make firts pnts[0] entry centre point

          x_c = reg.min_x + (i / step) -1;
          y_c = reg.min_y  + (i % step );
 

          // Pt ordering     	1 2 3
         //		4 0 8
         //		7 6 5

	  

     if (!get_point(&point,x_c , 0, y_c, 0, reg.subsampling, reg.image) )
       {  k1[i] = RIDICULOUS_NUMBER;
          continue;
       } */

      rnew->reset_all();
      ki = 0;
      for (ni = reg.minx; (ni <= reg.maxx) && ki < i; ni++)       
        for (nj = reg.miny; (nj <= reg.maxy) && ki < i; nj++)
        {  
          if (reg.included(ni,nj)) ki++;
          if (i == ki)
	  { ix = ni;
            iy = nj;
	    }     
	  }

      rnew->set_point(ix,iy);

      // Add three times neighborhood to obtain resonable
      // projection plane

       for (ki = 0; ki < 3; ki++)
        *rnew |= rnew->neighbourhood(); 

      // prevent *rnew from having points that have no data here

      *rnew &= reg; 
      
      pp1 = new pplane(*rnew,ind);
      delete pp1;        

      // assign indices to be of rnew to be the ones for reg
     
      for (j = 0; j < 9; j++) converted[j] = 0;
      
      ki = kj = 0;      

      for (ni = reg.minx; ni <= reg.maxx; ni++)
        for (nj = reg.miny; nj <= reg.maxy; nj++)
	{ if (reg.included(ni,nj)) ki++;
          if (rnew->included(ni,nj)) kj++;
          for (j = 0; j < 9; j++)
	    if (ind[j] == kj && !converted[j])
	    { ind[j] = ki;
              converted[j] = 1;
	      }
	  }

      index = ind[4];

      pnts[0].p[0] = list_x[index][1];
      pnts[0].p[1] = list_x[index][2];
      pnts[0].p[2] = list_y[index];

    
      // top left point

      /* if (!get_point(&point,x_c , -x_off , y_c , -y_off, reg.subsampling, reg.image) )
       {  k1[i] = RIDICULOUS_NUMBER;
          continue;
       } */

      index = ind[0];

      pnts[1].p[0] =  list_x[index][1];
      pnts[1].p[1] = list_x[index][2];
      pnts[1].p[2] =  list_y[index];

     // top middle point

      /* if (!get_point(&point, x_c , 0  , y_c, -y_off, reg.subsampling, reg.image) )
       {  k1[i] = RIDICULOUS_NUMBER;
          continue;
       } */
       
       index = ind[1];

       pnts[2].p[0] =  list_x[index][1];
       pnts[2].p[1] = list_x[index][2];
       pnts[2].p[2] = list_y[index];

      // top right point

       /* if (!get_point(&point, x_c , x_off , y_c,  y_off, reg.subsampling, reg.image) )
       {  k1[i] = RIDICULOUS_NUMBER;
          continue;
       } */

      index = ind[2];

      pnts[3].p[0] =  list_x[index][1];
      pnts[3].p[1] = list_x[index][2];
      pnts[3].p[2] =  list_y[index];

      // middle left point

      /* if (!get_point(&point, x_c,   -x_off , y_c, 0, reg.subsampling, reg.image) )
       {  k1[i] = RIDICULOUS_NUMBER;
          continue;
       } */

      index = ind[3]; 

      pnts[4].p[0] =  list_x[index][1];
      pnts[4].p[1] =  list_x[index][2];
      pnts[4].p[2] =  list_y[index];

      // middle right point
 
      /* if (!get_point(&point, x_c, x_off , y_c, 0, reg.subsampling, reg.image) )
       {  k1[i] = RIDICULOUS_NUMBER;
          continue;
       } */

      index = ind[5];

      pnts[8].p[0] =  list_x[index][1];
      pnts[8].p[1] =  list_x[index][2];
      pnts[8].p[2] =  list_y[index];

     // bottom left point
      /* if (!get_point(&point, x_c,   -x_off , y_c, y_off, reg.subsampling, reg.image) )
       {  k1[i] = RIDICULOUS_NUMBER;
          continue;
       } */

      index = ind[6];

      pnts[7].p[0] = list_x[index][1];
      pnts[7].p[1] = list_x[index][2];
      pnts[7].p[2] = list_y[index];

     // bottom middle point

      /* if (!get_point(&point, x_c, 0, y_c,  y_off, reg.subsampling, reg.image) )
       {  k1[i] = RIDICULOUS_NUMBER;
          continue;
       } */

      index = ind[7];

      pnts[6].p[0] =  list_x[index][1];
      pnts[6].p[1] = list_x[index][2];
      pnts[6].p[2] = list_y[index];


 // bottom right point

      /* if (!get_point(&point, x_c, x_off , y_c, y_off, reg.subsampling, reg.image) )
       {  k1[i] = RIDICULOUS_NUMBER;
          continue;
       } */

      index = ind[8];

      pnts[5].p[0] = list_x[index][1];
      pnts[5].p[1] = list_x[index][2];
      pnts[5].p[2] = list_y[index];

    


     // make n normal at centre point
    
    index = index2;

    normal[0] = norm[i][1];
    normal[1] = norm[i][2];
    normal[2] = norm[i][3];
     



      error=rcad_circurv(9,inds,pnts[0].p, normal, pnts, &k1[i], &k2[i], dir1, dir2);

     // printf("Circurv: %d k1: %lf k2: %lf\n", error, k1[i], k2[i]);

      // find largetst k1 or k2 and set this  as estimate for torus's k

     if (error < 0) 
     { k1[i] = k2[i] = RIDICULOUS_NUMBER;
       }
      
     k1[i] *= -1.0;
     k2[i] *= -1.0;

     if ( fabs(k1[i])  > fabs(k2[i]) )
      {  // sort so that k1 is always the smallest  
         temp = k1[i];
         k1[i] = k2[i];
         k2[i] = temp;
      }
    }

// now check curvature in k1 and k2 arrays

//  k1 should be constant, k2 can vary;   get the  pairs of k1, k2 to vote for torus type
sumk1 = sumk1sq = sumk2 = sumk2sq = 0.0;

  index = index_base;

torus_vote[0] = torus_vote[1] = 0;
step_no =0;

   for  (i=1;i<=no; i += step)
      {  if (k1[i] != RIDICULOUS_NUMBER)
         {  sumk1 += k1[i];
            sumk2 += k2[i];

            sumk1sq += k1[i]*k1[i];
            sumk2sq += k2[i]*k2[i];

           if   (  ( fabs(k1[i]) > fabs(k2[i]) ) ||   ( k1[i]*k2[i] < 0.0)  )  
              {  // apple torus suspected
                 ++torus_vote[1];
              }
          else  
              ++torus_vote[0];

           ++step_no;
        }
      }

     sumk1 /= step_no;
     sumk2 /= step_no;

     vark1 =  sumk1sq/step_no  - sumk1*sumk1;
     vark2 =  sumk2sq/step_no  - sumk2*sumk2;


     /* if (k1[index] == RIDICULOUS_NUMBER)
      {  // curvature not calculate for index pt
	//  so assign average value
        k1[index] = sumk1;
        k2[index] = sumk2;
      } */

     // added by Bojan Kverh 8.7.98

     // end added

     if ( vark1 > KTHRESH )
       { printf("Cannot Cope with this as a torus yet (KTHRESH),  k = %f,  s = %f, vark1 = %f, vark2 = %f\n",k1[index],k2[index],vark1,vark2);
          return (0);
       }
      else
       {  // assign torus params


          if (  fabs(k2[index]) < CYL_THRESH )               
             { // singularity needs addressing
               // above sorting will always place k2[index] as zero.

               // NEED to detect if patch is close to a cylinder as this need a special case
               // otherwise simply swap  k1 and k2 and use general torus fit as this is a cylinder and this      	                            
              // works ok

              // Check for cylinder   cross product  patch normals with axis estimate if below threshhold 
              //  (close to zero) the cylinder not cone.
               sumk1 = sumk1sq = 0.0;

             for (i = 1; i <= no;++i)
                 {    angle =   dot(norm[i],(atorus-1));
                      sumk1 += angle;
                      sumk1sq += angle*angle;
                 }
             sumk1 /= no;

              vark1 =  sumk1sq/no  - sumk1*sumk1;


	      if ( sumk1sq < CYL_THRESH )                    
               {  // cylinder
                   b[6] = k1[index];   // torus k
                   b[7] = k2[index];   // torus s
               }
           else
                { // cone type  (well at least not cylindrical)
                    
                      b[6] = k2[index];   // torus k
                      b[7] = k1[index];   // torus s
                 }

                   

                

             }
         else
          { // stay with above estimates
            // added by Bojan Kverh 8.7.98
            if (fabs(k1[index]) < fabs(k2[index]))
	    { b[7] = k1[index];
              b[6] = k2[index];
	      } else          // end added
            { b[6] = k1[index];   // torus k
	      b[7] = k2[index];   // torus s  
	      }
         }
          if ( torus_vote[0] > torus_vote[1] )
             b[8] = -1;  // lemon torus
           else
              b[8] = 1;  // apple torus

       }


    }
else

{  /* copy current fit parameters to b and q */

         b[1] = param[1];
         b[2] = param[2];
         b[3] = param[3];
         b[4] = param[4];
         b[5] = param[5];
         b[6] = param[6];
         b[7] = param[7];
         b[8] = param[8];

         q[1] = param[11];
         q[2] = param[12];
         q[3] = param[13];


       

         }

           for(i=1; i<=no; i++) {
                     /* set NR variables */
                    F[i] = 0.0;
                    sig[i] = 1.0;

             
              
             /* translate points to now coord system */

             list[i][1] = list_x[i][1] - q[1];
             list[i][2] = list_x[i][2] - q[2];
             list[i][3]=  list_y[i] - q[3];


        }
      
         
           lista = (int *) ivector(1,PARAM_SIZE);
           
           lista[1] = 1;
           lista[2] = 2;
           lista[3] = 3;
           lista[4] = 4;
           lista[5] = 5;
           lista[6] = 6;
           lista[7] = 7;
           lista[8] = 8;
           
           alamda = -1.0;
           newchisq = oldchisq = 1.0;

           printf("Initial Torus Estimates: %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f, %4.2f %4.2f\n", b[1],b[2],b[3],b[4],b[5], b[6], b[7], b[8]);

   // non-linear least squares estimate torus NOW at last

if ( fabs(b[7])  > CYL_THRESH) 
{  // standard torus fitting routines

         if (b[7] > 10000.0 || b[8] > 10000.0) printf("\n\n\n%lf TORUS SUCKS!\n\n\n\n\n",param[1]);

         o_mrqmin(list, F, sig, no,b,PARAM_SIZE, lista,PARAM_SIZE-1, covar, alpha, &newchisq,&alamda,  torusdx);
           
           k= 0;
           
           while (  k++ < 25 )
            { // iterate until no better error or max of 25 steps 
              oldchisq = newchisq;
              o_mrqmin(list, F, sig, no,b,PARAM_SIZE, lista,PARAM_SIZE-1, covar, alpha, &newchisq,&alamda,torusdx);
              
                           
            }
            
           alamda = 0.0;
           o_mrqmin(list, F, sig, no,b,PARAM_SIZE, lista,PARAM_SIZE-1, covar, alpha, &newchisq,&alamda,torusdx);
           
                         

	   // added by Bojan Kverh 1.6.98.
	   // modify the axis if necessary          
	    
  // if (b[6] * b[7] < 0.0) b[7] = - b[7]; 

         // end added by Bojan Kverh 1.6.1998
       
           param[1] = b[1];
           param[2] = b[2];
           param[3] = b[3]; 
           param[4] =  b[4] ; 
           param[5] = b[5];
           param[6] = b[6];
           param[7] = b[7];
           param[8] = b[8];

           param[11]  = q[1];
           param[12] = q[2];
           param[13] = q[3];


          }
else
{  // degenerat case fit cylinder for now but allow patch reestimates by setting b[8] (epsilon) to 0
  //  this is checked for above for first estimate call

        printf("Torus Fitting: Degenerate case use Cylinder fit for now\n");

        b[8] = 0;

        b[1] = RIDICULOUS_NUMBER;

          
          

          // fit away

          cyl_fit(list_x, list_y, norm, no, b, chisq);

          // map back to suit torus way of params;

          
           param[1] = b[1];
           param[2] = b[2];
           param[3] = b[3]; 

           
           /* sigma = atan2(atorus[2]/atorus[1]);  but in NORAML C ARRAY INDEX!! */

       param[4] = b[4];   // this is now alpha


     /* tau = acos(n[3]) but in NORMAL C ARREAY INDEX!! */

       param[5] =  0.0; // NOT used

           
           param[6] = b[5];   //  k = cyl k
           param[7] = b[7];
           param[8] = b[8];

           param[11]  = q[1];
           param[12] = q[2];
           param[13] = q[3];

      }   

        

     

     // tidy up nr memory allocs


    free_dvector(b,1,PARAM_SIZE);
    free_dvector(n,1,3);
    free_dmatrix(list,1,no,1,3);
    free_dmatrix(covar,1,PARAM_SIZE,1,PARAM_SIZE);
    free_dmatrix(alpha,1,PARAM_SIZE,1,PARAM_SIZE);
  
    free_ivector(lista,1,PARAM_SIZE);
    //    free_dvector(b,1,PARAM_SIZE);
    free_dvector(F,1,no);
    free_dvector(sig,1,no );

    free_dvector(k1,1,no);

    free_dvector(k2,1,no);


    free_dvector(q,1,3);

return(1);
}

/*
int torus::get_point(struct vect *point, int x, int x_off, int y, int y_off, int subsampling, RangeImage  *img)


{     int im_x = x*subsampling + x_off, im_y = y*subsampling + y_off;

      if (im_x >= 0 && im_x < img->rwidth() && im_y >= 0 && im_y < img->rheight())
          {   if ( (point->z =  img->imz(im_x,im_y)) != 0.0 )
                { point->x =  img->imx(im_x,im_y);
                  point->y =  img->imy(im_x,im_y);
   
                  return(1);
               }
             else
               return(0); // no depth value
          }
    else 

      return(0);
}
*/


void torus::fprint(FILE *f)
{ fprintf(f,"7 %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],
    a[11],a[12],a[13]);
  }


void torus::parameters(char **name, double *value)
{ int i;
  sprintf(name[0],"rho");
  sprintf(name[1],"phi");
  sprintf(name[2],"theta");
  sprintf(name[3],"sigma");
  sprintf(name[4],"tau");
  sprintf(name[5],"k");
  sprintf(name[6],"s");
  sprintf(name[7],"epsilon");
  for (i = 0; i < 8; i++)
    value[i] = a[i+1];

  value[8] = a[11];
  value[9] = a[12];
  value[10] = a[13];
}

void torus::set_parameters(double *value)
{ int i;
  for (i = 0; i < 8; i++)
    a[i+1] = value[i];

  a[11] = value[8];
  a[12] = value[9];
  a[13] = value[10];
  }
  
struct point torus::project(struct point &p)
{ vector3D axis,n,axis2,n2,x,v,axis3,r2,r,c,mr;
  double kot,vr,h,nrm;  
  struct point p1;

  axis.el(0) = cos(a[4])*sin(a[5]);
  axis.el(1) = sin(a[4])*sin(a[5]);
  axis.el(2) = cos(a[5]);

  n.el(0) = cos(a[2])*sin(a[3]);
  n.el(1) = sin(a[2])*sin(a[3]);
  n.el(2) = cos(a[3]);

  v.el(0) = p.x - a[11];
  v.el(1) = p.y - a[12];
  v.el(2) = p.z - a[13];

  kot = acos(n.inner(axis));

  vr = fabs(1/a[7]-1/a[6])*sin(kot);
  h = fabs(1/a[7]-1/a[6])*cos(kot);
  
  axis2 = axis;
  axis2.multiply_col(0,h);
  n2 = n;
  n2.multiply_col(0,a[1]+1/a[7]);

  x = v - n2 - axis2;

  std::cout << "Axis: " << axis << "  N: " << n;
  std::cout << "VR = " << vr << "  height: " << h;

  c = n2 + axis2;
  c.el(0) += a[11];
  c.el(1) += a[12];
  c.el(2) += a[13];

  std::cout << "Center of torus " << c;

  axis3 = axis.outer(x);
  axis3 = axis3.normalize();
  r2 = axis3.outer(axis);
  r2 = r2.normalize();

  r2.multiply_col(0,vr);
  mr = v - n2 - axis2;
  std::cout << "x = " << mr << "  R2 = " << r2 << " Relation: " << mr.el(0) / mr.el(2) << "  " << r2.el(0)/r2.el(2);
  r = v - n2 - axis2 - r2;
  nrm = r.norm();
  
  r = r.normalize();
  r.multiply_col(0,nrm-1/a[6]);

  p1.x = p.x - r.el(0);
  p1.y = p.y - r.el(1);
  p1.z = p.z - r.el(2);
  
  std::cout << "\nProjected point distance: " << distance(p1);
  std::cout << "\nNormal point distance: " << distance(p) << "\n\n";

  return(p1);
  }

