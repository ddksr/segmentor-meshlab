// cone.C

#include "image.h"
#include "matrix.h"
#include "cone.h"

/* Ralph`c curvature fitting code */

#include "pplane.h"


#include <math.h>
#include <stdlib.h>



extern "C" {

  #include "rcad_vector.h"
  #include "rcad_ndm_ut.h"
  #include "rcad_circurv.h"
  
void o_mrqmin(double **x, double *y,double *sig,int ndata,double *a,int ma,int *lista,int mfit,double **covar,double **alpha,double *chisq,double *alamda, int (*funcs)(double,double ,double, double *,double *,double *,int));

}

#include "rcad_circurv.h"
  
// extern "C" int rcad_circurv(int, int *, double *, double *, rcad_SurfaceData *, double *, double *, double *, double *);

#define TRUE 1
#define FALSE 0

#define PARAM_SIZE 6


#define RIDICULOUS_NUMBER -1000000000

// double cone::m_error = max_description_error;   // provide maximum point distance and description error
// double cone::m_dist = max_point_distance;    // for all objects of the class

extern double gmax_description_error, gmax_point_distance;


cone::cone(region &reg) {
 
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
  
  // m_error = gmax_description_error;   // provide maximum point distance and description error
  // m_dist = gmax_point_distance;    // 

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
  	  std::cout << "Over the limit of Cylinder points";

  param[1] = RIDICULOUS_NUMBER;

int xspan = max_x - min_x;
int yspan = max_y - min_y;

if( cone_fit(reg, list_x, list_y, list_norm, no, xspan, yspan, param, &chisq) )
 
   { // Gabor rep 
     a[1] = param[1];
     a[2] = param[2];
     a[3] = param[3];
     a[4] = param[4];
     a[5] = param[5];
     a[6] = param[6];

     a[11] = param[11];
     a[12] = param[12];
     a[13] = param[13];
     adjust_z(reg);
}

 else
 { // make params force an large error

    a[1] = RIDICULOUS_NUMBER;
 }
 
 
  fflush(stdout);

  free_dmatrix(list_x,1,Max,1,2);
  free_dmatrix(list_norm,1,Max,1,3);

  free_dvector(list_y,1,Max);
  free_dvector(param,1,13);

}

cone::cone(region &reg, cone &init) {
 
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
  
  // m_error = gmax_description_error;   // provide maximum point distance and description error
  // m_dist = gmax_point_distance;    // 

     param[1] = init.a[1];
     param[2] = init.a[2];
     param[3] = init.a[3];
     param[4] = init.a[4];
     param[5] = init.a[5];
     param[6] = init.a[6];


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

if( cone_fit(reg, list_x, list_y, list_norm, no, xspan, yspan, param, &chisq) )
 
   { /* Gabor rep */
     a[1] = param[1];
     a[2] = param[2];
     a[3] = param[3];
     a[4] = param[4];
     a[5] = param[5];
     a[6] = param[6];

     a[11] = param[11];
     a[12] = param[12];
     a[13] = param[13];
     adjust_z(reg); 
}

else
 { // make params force an large error

    a[1] = 1000000000.0;
 }

     
 
  free_dmatrix(list_x,1,Max,1,2);
  free_dvector(list_y,1,Max);
  free_dvector(param,1,13);
  free_dmatrix(list_norm,1,Max,1,3);

}

// Euclidian distance function - Implemented by Bojan Kverh 21.5.98

double cone::abs_signed_distance(struct point &v)
{ vector3D axis,n,acn,c,axis3,n2,origin;
  double pcd,pca,d,kot;

  axis.el(0) = cos(a[4]) * sin(a[5]);
  axis.el(1) = sin(a[4]) * sin(a[5]);
  axis.el(2) = cos(a[5]);
  n.el(0) = cos(a[2]) * sin(a[3]);
  n.el(1) = sin(a[2]) * sin(a[3]);
  n.el(2) = cos(a[3]);

  origin.el(0) = a[11];
  origin.el(1) = a[12];
  origin.el(2) = a[13];

  axis3 = axis;
  if (a[6] > 0.0) axis3.multiply_col(0,-1);
  kot = acos(n.inner(axis3));
  axis3 = axis;
  axis3.multiply_col(0,(1/fabs(a[6]))/cos(kot));
  n2 = n;
  n2.multiply_col(0,a[1]+1/a[6]);
  c = origin + n2 + axis3;
  c.el(0) = v.x - c.el(0);
  c.el(1) = v.y - c.el(1);
  c.el(2) = v.z - c.el(2);

  pcd = c.norm();
  pca = c.inner(axis);

  acn = n.outer(axis);

  d = acn.norm()*sqrt(pcd*pcd-pca*pca)-fabs(n.inner(axis)*c.inner(axis));

  return(d);
  }

/*
double cone::distance(struct point& v)
{  double r_est;
   double modsqncrossa;  

    double sinphi = sin(a[2]), cosphi = cos(a[2]);
    double sintheta = sin(a[3]), costheta = cos(a[3]);
    double sinsigma = sin(a[4]), cossigma = cos(a[4]);
    double sintau = sin(a[5]), costau = cos(a[5]);

      
    //   double *p  = dvector(1,3);
    //   double *phat  = dvector(1,3);
    //   double *acone  = dvector(1,3);
    //   double *n  = dvector(1,3);
    //   double *ncrossa = dvector(1,3);

   vector3D p,phat,acone,n,ncrossa;

   double phatdota, phatdotn, ndota;

    double top,bottom;


    p.el(0) = v.x - a[11];
    p.el(1) = v.y - a[12];
    p.el(2) = v.z - a[13];

    // calculate n

    n.el(0) = cosphi*sintheta;
    n.el(1) = sinphi*sintheta;
    n.el(2) = costheta;

    // phat - p - rho*n

     phat.el(0) = p.el(0) - a[1]*n.el(0);
     phat.el(1) = p.el(1) - a[1]*n.el(1);
     phat.el(2) = p.el(2) - a[1]*n.el(2);

    // compute acone

     acone.el(0) = cossigma*sintau;
     acone.el(1) = sinsigma*sintau;
     acone.el(2) = costau;

     // compute dot products

     phatdota = phat.inner(acone);
     phatdotn = phat.inner(n);
     ndota = n.inner(acone);

     // compute n cross a product

    ncrossa = n.outer(acone);
    modsqncrossa = ncrossa.norm();

    top = a[6]*(modsqncrossa*phat.norm() - phatdota*phatdota)*0.5 - phatdotn*modsqncrossa;

    bottom = a[6]*phatdota*ndota + modsqncrossa;

    if (fabs(bottom) > 0.0000000001) r_est = top/bottom;
      else r_est = RIDICULOUS_NUMBER; 
    // tidy up memory of nr utils

    //    free_dvector(ncrossa,1,3);
    //    free_dvector(p,1,3);
    //    free_dvector(phat,1,3);
    //    free_dvector(n,1,3);
    //    free_dvector(acone,1,3);




    return fabs(r_est);


        
 }
*/

void cone::save(FILE *f) { 
  fwrite(a,14,sizeof(double),f);
  fwrite(&zmin,1,sizeof(double),f);
  fwrite(&zmax,1,sizeof(double),f);
  
}

cone::cone(FILE *f) {
  fread(a,14,sizeof(double),f);
  fread(&zmin,1,sizeof(double),f);
  fread(&zmax,1,sizeof(double),f);
 
}

model* cone::improve(region& orig, region& added)
{ region r = orig | added;
  return(new cone(r));
  }

model* cone::improve(model *m, region& orig, region& added)
{ region r = orig | added;
  return(new cone(r,*(cone *)m));
  }


int conedx(double x, double y, double z, double param[], double *f, double dfda[], int na)

{  double det;
     
    double sinphi = sin(param[2]), cosphi = cos(param[2]);
    double sintheta = sin(param[3]), costheta = cos(param[3]);
    double sinsigma = sin(param[4]), cossigma = cos(param[4]);
    double sintau = sin(param[5]), costau = cos(param[5]);

   double modsqphat,modsqncrossa;  
      
   double *p  = dvector(1,3);
   double *p2  = dvector(1,3);
   double *phat  = dvector(1,3);
   double *acone  = dvector(1,3);
   double *dadsigma  = dvector(1,3);
   double *dadtau  = dvector(1,3);
   double *n  = dvector(1,3);
   double *dndphi = dvector(1,3);
   double *dndtheta  = dvector(1,3);
   double *ncrossa = dvector(1,3);

   double pdota, pdotn,
              phatdota, phatdotn, 
              phatdotdadsigma, phatdotdadtau,
              phatdotdndphi, phatdotdndtheta,
              ndota,  ndotdadsigma, ndotdadtau,
              dndphidotp,dndphidota,
              dndthetadotp,dndthetadota;
   
   double p2dota;

   double lambda,xi, mu, eta;

    double dlambdadrho, dlambdadphi, dlambdadtheta, dlambdadsigma, dlambdadtau;
    double dxidrho, dxidphi, dxidtheta, dxidsigma, dxidtau;
    double dmudrho, dmudphi, dmudtheta, dmudsigma, dmudtau;
    double detadrho, detadphi, detadtheta, detadsigma, detadtau;

    double top,bottom;

  
    if (na) na = na;  // to prevent stupid compiler warnings

    p[1] = x;
    p[2] = y;
    p[3] = z;

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

    // p2 = p - 2*rho*n

     p2[1] = p[1] - 2.0*param[1]*n[1];
     p2[2] = p[2] - 2.0*param[1]*n[2];
     p2[3] = p[3] - 2.0*param[1]*n[3];

    // compute acone, dadsigma, dadtau

     acone[1] = cossigma*sintau;
     acone[2] = sinsigma*sintau;
     acone[3] = costau;

     dadsigma[1] = -sinsigma*sintau;
     dadsigma[2] = cossigma*sintau;
     dadsigma[3] = 0;

     dadtau[1] = cossigma*costau;
     dadtau[2] = sinsigma*costau;
     dadtau[3] = -sintau;


     // compute dot products
     pdota = dot(p,acone);
     pdotn = dot(p, n);

     phatdota = dot(phat,acone);
     phatdotn = dot(phat, n);
     phatdotdadsigma = dot(phat, dadsigma);
     phatdotdadtau = dot(phat, dadtau);
     phatdotdndphi = dot(phat, dndphi);
     phatdotdndtheta = dot(phat, dndtheta);

     p2dota = dot(p2,acone);

     ndota = dot(n,acone);
     ndotdadsigma =  dot(n,dadsigma);
     ndotdadtau =  dot(n,dadtau);

     dndphidotp = dot(dndphi,p);
     dndphidota = dot(dndphi,acone);

     dndthetadotp = dot(dndtheta,p);
     dndthetadota = dot(dndtheta,acone);


     // compute n cross a product

    cross_prod(n,acone,ncrossa);
    modsqncrossa = modsq(ncrossa);

    // compute lambda,xi, mu, eta and their derivatives
 
    lambda = (modsqncrossa*modsqphat - phatdota*phatdota)*0.5;
    xi = - phatdotn*modsqncrossa;
    mu = phatdota*ndota;
    eta = modsqncrossa;

    dlambdadrho = param[1]*(modsqncrossa - 2*ndota*ndota) + 2*pdota*ndota - pdotn*modsqncrossa;
    dlambdadphi =  -1.0*dndphidota*ndota*modsqphat - param[1]*modsqncrossa*dndphidotp + 2*param [1]*phatdota*dndphidota;
    dlambdadtheta = -1.0*dndthetadota*ndota*modsqphat - param[1]*modsqncrossa*dndthetadotp + 2*param[1]*phatdota*dndthetadota;
      dlambdadsigma = -1.0*(ndotdadsigma*ndota*modsqphat + phatdota*phatdotdadsigma);
     dlambdadtau = -1.0*(ndotdadtau*ndota*modsqphat + phatdota*phatdotdadtau);

     dxidrho = 2.0*modsqncrossa;
     dxidphi = 2.0*dndphidota*ndota - phatdotdndphi*modsqncrossa;
     dxidtheta = 2.0*dndthetadota*ndota - phatdotdndtheta*modsqncrossa;
     dxidsigma = phatdotn*ndotdadsigma*ndota;
     dxidtau = phatdotn*ndotdadtau*ndota;

     dmudrho = -1.0*modsqncrossa;
     dmudphi =  dndphidota*p2dota;
     dmudtheta = dndthetadota*p2dota;
     dmudsigma = phatdotdadsigma*ndota + phatdota*ndotdadsigma;
     dmudtau = phatdotdadtau*ndota + phatdota*ndotdadtau;

     detadrho = 0.0;
     detadphi = -2.0*dndphidota*ndota;
     detadtheta = -2.0*dndthetadota*ndota;
     detadsigma = -2.0*ndotdadsigma*ndota;
     detadtau = -2.0*ndotdadtau*ndota;

     det = (mu*param[6] + eta);
     det = det*det;

      /*   compute fit with current params */

      top = param[6]*(modsqncrossa*modsq(phat) - phatdota*phatdota)*0.5 - phatdotn*modsqncrossa;

      bottom = param[6]*phatdota*ndota + modsqncrossa;

      *f = top/bottom;

      /* d/drho */   
     dfda[1] = ( (dlambdadrho*mu - lambda*dmudrho)*param[6]*param[6] +
                      (dlambdadrho*eta + mu*dxidrho - xi*dmudrho)*param[6] + dxidrho*eta) / det;

      /* d/dphi */

     dfda[2] =   ( (dlambdadphi*mu - lambda*dmudphi)*param[6]*param[6] +
                       (dlambdadphi*eta + mu*dxidphi - lambda*detadphi - xi*dmudphi)*param[6]
                       + dxidphi*eta - xi*detadphi) / det;

     /* d/dtheta */

     dfda[3] =    ( (dlambdadtheta*mu - lambda*dmudtheta)*param[6]*param[6] +
                         (dlambdadtheta*eta + mu*dxidtheta - lambda*detadtheta - xi*dmudtheta)*param[6]
                       + dxidtheta*eta - xi*detadtheta) / det;
  
      /* d/dsigma */

     dfda[4] =   ( (dlambdadsigma*mu - lambda*dmudsigma)*param[6]*param[6] +
                         (dlambdadsigma*eta + mu*dxidsigma - lambda*detadsigma - xi*dmudsigma)*param[6]
                       + dxidsigma*eta - xi*detadsigma) / det;

     /* d/dtau */

     dfda[5] =   ( (dlambdadtau*mu - lambda*dmudtau)*param[6]*param[6] +
                         (dlambdadtau*eta + mu*dxidtau - lambda*detadtau - xi*dmudtau)*param[6]
                       + dxidtau*eta - xi*detadtau) / det;

    /* d/dk */

     dfda[6] =    (lambda*eta - mu*xi) / det;
  
   // tidy up memory of nr utils

    free_dvector(p,1,3);
    free_dvector(p2,1,3);
    free_dvector(phat,1,3);
    free_dvector(acone,1,3);
    free_dvector(dadsigma,1,3);
    free_dvector(dadtau,1,3);
    free_dvector(n,1,3);
    free_dvector(dndphi,1,3);
    free_dvector(dndtheta,1,3);
    free_dvector(ncrossa,1,3);

    return 1;
 
}

int cone::cone_fit(region &reg, double **list_x, double *list_y, double **norm, int no, int xspan, int yspan,double *param, double *chisq)

{
  double pt[3], pt_norm[3];
  double  acone[3];

  int ind[9];         // indices of  extreme points

  sline input_norms[9];
  sline output_axis;
  sline *r_axis=&output_axis;
  int   axnum=1;

  int i, index, k;

   if (chisq) chisq = chisq;    // to prevent stupid compiler warnings.
   
   double *q= (double *) dvector(1,3);  /* cone coord sys origin */

 
 // Numerical recipe vars
     
     double **list, *F, *sig;
     double  *b = dvector(1,PARAM_SIZE);
     double  *n = dvector(1,3);
      int *lista;
      double alamda, newchisq, oldchisq;

      
      double  **covar = (double **) dmatrix(1,PARAM_SIZE,1,PARAM_SIZE),
                  **alpha = (double **) dmatrix(1,PARAM_SIZE,1,PARAM_SIZE); 


      list = (double **) dmatrix(1,no,1,3);
      F = (double *) dvector(1,no);
      sig = (double *) dvector(1,no);


      pplane projection_plane(reg,ind);           // fills ind with indices of extremes

if ( param[1] == RIDICULOUS_NUMBER )
  { /* NEW REGION SO CONSTRUCT INITIAL ESTIMATES */

    
  // Guess AXIS of CONE 

  //  Take 5 pts at 4 corners and center of region and compute sline of normals

    // top left normal

    pt[0] =  list_x[ind[0]][1];
    pt[1] = list_x[ind[0]][2];
    pt[2] = list_y[ind[0]];

    for (i=0;i<3;++i)
      pt_norm[i] = norm[ind[0]][i+1];

   input_norms[0] = sline(pt, pt_norm,0); 

    // top right normal

    // changed by Bojan Kverh 25.2.1998 

   index = ind[2];    // xspan+1;

    pt[0] =  list_x[index][1];
    pt[1] = list_x[index][2];
    pt[2] = list_y[index];

    for (i=0;i<3;++i)
      pt_norm[i] = norm[index][i+1];


   input_norms[1] = sline(pt, pt_norm,0);
 
   // middle normal

   index = ind[4];

    pt[0] =  list_x[index][1];
    pt[1] = list_x[index][2];
    pt[2] = list_y[index];

    for (i=0;i<3;++i)
      pt_norm[i] = norm[index][i+1];

   input_norms[2] = sline(pt, pt_norm,0); 

    // bottom left normal

    // changed by Bojan Kverh 25.2.1998 

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
    pt[2 ] = list_y[index];

    for (i=0;i<3;++i)
      pt_norm[i] = norm[index][i+1];

   input_norms[4] = sline(pt, pt_norm,0); 

   index = guess_axes(input_norms,5,
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

    output_axis.direction(acone);

 // index = 1;
 // acone[0] = acone[2] = 0.0;
//  acone[1] = 1.0;

 if ( (index == 0) || ( fabs(modsq(acone-1) - 0.999999) > 0.001 ))
       {  // we have no axis
         std::cout << "ERROR: Cannot Guess Axis for this Cone\n" << "index = " << index << "\n";
         std::cout << "acone[0] =" << acone[0]  << " acone[1] =" << acone[1]  << " acone[2] =" << acone[2] << "\n";
         std::cout << "modsq = " << modsq(acone-1) << "\n";
         return(0);

       }

  // printf("%d:%f %f %f\n",index,acone[0],acone[1],acone[2]);
  

  // compute  initial estimates

   // make n normal at first point

    n[1] = -norm[1][1];
    n[2] = -norm[1][2];
    n[3] = -norm[1][3];

 // make rho = 1

   b[1] = 1.0;

   // Make origin of cone coordinate system a point 1 "unit" off 
   //      1st cone point along normal at the point

      q[1] = list_x[1][1] -  b[1]* n[1];
      q[2] = list_x[1][2] -  b[1]* n[2];
      q[3] = list_y[1] -  b[1]* n[3];

      /*          
      // Added by Bojan Kverh 26.2.1998    
      // Estimation of cone axis
      // As this will be done only at the beginning,
      // the axis is very similar to the one of the cylinder
      
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
      
      acone[0] = ax.el(0) / npts;
      acone[1] = ax.el(1) / npts;
      acone[2] = ax.el(2) / npts;

      // End Added by Bojan Kverh 26.2.1998 
      */
      
   
      // compute estimates of angles

      // phi = atan2(n[2]/n[1]);

       b[2] = atan2(n[2],n[1]);


       // theta = acos(n[3])

         b[3] = acos(n[3]);


	 // sigma = atan2(acone[2]/acone[1]);  but in NORAML C ARREAY INDEX!!

       b[4] = atan2(acone[1],acone[0]);


       // tau = acos(n[3]) but in NORAML C ARREAY INDEX!!

       b[5] = acos(acone[2]);

       
     //  To Estimate k  --- use ralph`s curvature estimation.
    //    SInce COne one curvature should be zero or close to it.
     //   Choose no zero curvature value as estimate.
     
      rcad_SurfaceData pnts[9];
      int              inds[9];

      // make firts pnts[0] entry centre point

      // Pt ordering     	1 2 3
      //		4 0 8
      //		7 6 5

      xspan = yspan = (int)sqrt(no);
     
      index =  ind[4];
      
      pnts[0].p[0] = list_x[index][1];
      pnts[0].p[1] = list_x[index][2];
      pnts[0].p[2] = list_y[index];

    
      // top left point


      index = ind[0];

      pnts[1].p[0] =  list_x[index][1];
      pnts[1].p[1] = list_x[index][2];
      pnts[1].p[2] = list_y[index];

     // top middle point


      index = ind[1];

       pnts[2].p[0] =  list_x[index][1];
       pnts[2].p[1] = list_x[index][2];
       pnts[2].p[2] = list_y[index];

      // top right point


      index = ind[2];

      pnts[3].p[0] =  list_x[index][1];
      pnts[3].p[1] = list_x[index][2];
      pnts[3].p[2] = list_y[index];

      // middle left point


      index = ind[3];

      pnts[4].p[0] =  list_x[index][1];
      pnts[4].p[1] = list_x[index][2];
      pnts[4].p[2] = list_y[index];

      // middle right point


      index =  ind[5];

      pnts[8].p[0] =  list_x[index][1];
      pnts[8].p[1] = list_x[index][2];
      pnts[8].p[2] = list_y[index];

     // bottom left point

      index = ind[6];

      pnts[7].p[0] =  list_x[index][1];
      pnts[7].p[1] = list_x[index][2];
      pnts[7].p[2] = list_y[index];

     // bottom middle point

      index = ind[7];

      pnts[6].p[0] =  list_x[index][1];
      pnts[6].p[1] = list_x[index][2];
      pnts[6].p[2] = list_y[index];


 // bottom right point

    index = ind[8];

      pnts[5].p[0] =  list_x[index][1];
      pnts[5].p[1] = list_x[index][2];
      pnts[5].p[2] = list_y[index];

     for (i=0;i<9;++i)
       inds[i] = i;


     // make n normal at centre point
      index = ind[4];

     double normal[3];
     double dir1[3];
     double dir2[3];


    normal[0] = norm[index][1];
    normal[1] = norm[index][2];
    normal[2] = norm[index][3];
 

     double k1,k2;
     int error;

     error=rcad_circurv(9,inds,pnts[0].p, normal, pnts, &k1, &k2, dir1, dir2);

     //std::printf("Circurv: %d k1: %lf k2: %lf\n", error, k1, k2);

      // find largest k1 or k2 and set this  as estimate for cone's k


     if ( fabs(k1) > fabs(k2) )
       b[6] = -k1;
    else
       b[6] = -k2; 

}
else

{  /* copy current fit parameters to b and q */

         b[1] = param[1];
         b[2] = param[2];
         b[3] = param[3];
         b[4] = param[4];
         b[5] = param[5];
         b[6] = param[6];

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
           
           alamda = -1.0;
           newchisq = oldchisq = 1.0;

           //printf("Initial Cone Estimates: %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f\n", b[1],b[2],b[3],b[4],b[5], b[6]);

   // non-linear least squares estimate cone NOW at last


o_mrqmin(list, F, sig, no,b,PARAM_SIZE, lista,PARAM_SIZE, covar, alpha, &newchisq,&alamda,  conedx);
           
           k= 0;
           
           while (  k++ < 25 )
            { /* iterate until no better error or max of 25 steps */
              oldchisq = newchisq;
              o_mrqmin(list, F, sig, no,b,PARAM_SIZE, lista,PARAM_SIZE, covar, alpha, &newchisq,&alamda,conedx);
              
                           
            }
            
           alamda = 0.0;
           o_mrqmin(list, F, sig, no,b,PARAM_SIZE, lista,PARAM_SIZE, covar, alpha, &newchisq,&alamda,conedx);
           
	   // added by Bojan Kverh 21.5.98.
	   // modify the axis if necessary          
	   
 vector3D axis, n2, flag;
 double pi = 3.1415926535; 

         axis.el(0) = cos(b[4])*sin(b[5]);
         axis.el(1) = sin(b[4])*sin(b[5]);
         axis.el(2) = cos(b[5]);
         n2.el(0) = cos(b[2])*sin(b[3]);
         n2.el(1) = sin(b[2])*sin(b[3]);
         n2.el(2) = cos(b[3]);
         
         if (b[6] * n2.inner(axis) > 0.0) 
         { /* if (b[4] > 0) b[4] -= pi;
             else b[4] += pi; */
           if (b[5] > 0) b[5] -= pi;
             else b[5] += pi;
           }         
                 
           param[1] = b[1];
           param[2] = b[2];
           param[3] = b[3]; 
           param[4] =  b[4]; 
           param[5] = b[5];
           param[6] = b[6];

           param[11] = q[1];
           param[12] = q[2];
           param[13] = q[3];

     



     // tidy up nr memory allocs


    free_dvector(b,1,PARAM_SIZE);
    free_dvector(n,1,3);
    free_dmatrix(list,1,no,1,3);
    free_dmatrix(covar,1,PARAM_SIZE,1,PARAM_SIZE);
    free_dmatrix(alpha,1,PARAM_SIZE,1,PARAM_SIZE);
  
    free_ivector(lista,1,PARAM_SIZE);
    free_dvector(F,1,no);
    free_dvector(sig,1,no );

    free_dvector(q,1,3);

return(1);
}

void cone::adjust_z(region &reg)
{ vector3D c,p,axis,n,axis2,n2,origin;
  int i,j;
  double beta,d;
  struct point pt;

  axis.el(0) = cos(a[4]) * sin(a[5]);
  axis.el(1) = sin(a[4]) * sin(a[5]);
  axis.el(2) = cos(a[5]);
  n.el(0) = cos(a[2]) * sin(a[3]);
  n.el(1) = sin(a[2]) * sin(a[3]);
  n.el(2) = cos(a[3]);
  
  origin.el(0) = a[11];
  origin.el(1) = a[12];
  origin.el(2) = a[13];

  // cout << "\n" << a[4] << "  " << a[5] << "\n";

  // cout << axis;

  axis2 = axis;
  if (a[6] > 0.0) axis2.multiply_col(0,-1);
  beta = acos(n.inner(axis2));
  axis2 = axis;
  axis2.multiply_col(0,(1/fabs(a[6]))/cos(beta));
  n2 = n;
  n2.multiply_col(0,a[1]+1/a[6]);
  c = origin + n2 + axis2;

  zmin = 1e40;
  zmax = -zmin;

  // cout << axis;

  for (i = reg.minx; i <= reg.maxx; i++)
    for (j = reg.miny; j <= reg.maxy; j++)
      if (reg.included(i,j))
      { pt = reg.get_point(i,j);
        p.el(0) = pt.x - c.el(0);
        p.el(1) = pt.y - c.el(1);
        p.el(2) = pt.z - c.el(2);
        d = p.inner(axis);
        if (d < zmin) zmin = d;
        if (d > zmax) zmax = d;
	}

  if (zmin > 0 || zmax > 0)
  { if (zmin > 0 && zmax < 0 || zmin < 0 && zmax > 0)
      if (fabs(zmin) > fabs(zmax)) zmax = 0;
        else zmin = 0;
    if (zmin > 1e-6 || zmax > 1e-6)
    { a[4] = a[4];
      a[5] = a[5];
      d = zmax;
      zmax = -zmin;
      zmin = -d;
      }
    }
  }

void cone::fprint(FILE *f)
{ fprintf(f,"6 %lf %lf %lf %lf %lf %lf %lf %lf %lf",a[1],a[2],a[3],a[4],a[5],a[6],a[11],a[12],a[13]);
  }

void cone::draw() 
{ switch(draw_type)
  { case M_XWINDOWS:
      break;
    case M_OPENGL:
      drawGL();
      break;
    default:
      break;
    }
  }

void cone::parameters(char **name, double *value)
{ int i;
  sprintf(name[0],"rho");
  sprintf(name[1],"phi");
  sprintf(name[2],"theta");
  sprintf(name[3],"sigma");
  sprintf(name[4],"tau");
  sprintf(name[5],"k");
  for (i = 0; i < 6; i++)
    value[i] = a[i+1];

  value[6] = a[11];
  value[7] = a[12];
  value[8] = a[13];
  value[9] = zmin;
  value[10] = zmax;

}

void cone::set_parameters(double *value)
{ int i;
  for (i = 0; i < 6; i++)
    a[i+1] = value[i];

  a[11] = value[6];
  a[12] = value[7]; 
  a[13] = value[8];
  zmin = value[9];
  zmax = value[10];
  }

void cone::rif_write(FILE *f)
{ vector3D axis,n,n2,n2crossa,add,b;
  double d = a[1] + 1 / a[6],cospsi,sinpsi,rad1,rad2; 

  axis.el(0) = cos(a[4])*sin(a[5]);
  axis.el(1) = sin(a[4])*sin(a[5]);
  axis.el(2) = cos(a[5]);
    
  // cout << "\nDraw axis: " << axis;

  n.el(0) = cos(a[2])*sin(a[3]);
  n.el(1) = sin(a[2])*sin(a[3]);
  n.el(2) = cos(a[3]);
  n2 = n;
  n.multiply_col(0,d);

  n2crossa = n2.outer(axis);
  
  cospsi = n2crossa.norm();

  add.el(0) = n.el(0) + a[11];
  add.el(1) = n.el(1) + a[12];
  add.el(2) = n.el(2) + a[13];

  b = add + axis;

  sinpsi = n2.inner(axis);

  rad1 = (1/a[6])/cospsi;
  rad2 = fabs(rad1 + sinpsi / cospsi);

  fprintf(f,"G 5 {%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf};\n",add.el(0),add.el(1),add.el(2),fabs(rad1),b.el(0),b.el(1),b.el(2),rad2);
  }

struct point cone::project(struct point &p)
{ vector3D axis,n,v,r,x;
  double d = a[1] + 1 / a[6], z1, kot,t,nrm; 
  struct point p1;

  axis.el(0) = cos(a[4])*sin(a[5]);
  axis.el(1) = sin(a[4])*sin(a[5]);
  axis.el(2) = cos(a[5]);
  n.el(0) = cos(a[2])*sin(a[3]);
  n.el(1) = sin(a[2])*sin(a[3]);
  n.el(2) = cos(a[3]);
  //if (a[6] > 0.0) axis.multiply_col(0,-1.0);
  if (a[6] * axis.inner(n) > 0.0 ) axis.multiply_col(0,-1.0);

  kot = acos(-n.inner(axis));
  n.multiply_col(0,d);

  v.el(0) = p.x - a[11];
  v.el(1) = p.y - a[12];
  v.el(2) = p.z - a[13];

  x = v - n;
  z1 = axis.inner(x);  
  t = ((1/a[6])/sin(kot)) - z1 / tan(kot);

  axis.multiply_col(0,z1);
  r = v - (n + axis);
  nrm = r.norm();
  r = r.normalize();
  r.multiply_col(0,t - nrm);

  p1.x = p.x + r.el(0);
  p1.y = p.y + r.el(1);
  p1.z = p.z + r.el(2);

  return(p1);
  }

void cone::rotate(symatrix& A)
{ class vector axis(4),n(4),o(4);

  axis.el(0) = cos(a[4])*sin(a[5]);
  axis.el(1) = sin(a[4])*sin(a[5]);
  axis.el(2) = cos(a[5]);
  n.el(0) = cos(a[2])*sin(a[3]);
  n.el(1) = sin(a[2])*sin(a[3]);
  n.el(2) = cos(a[3]);
  n.el(3) = 1.0;
  // if (a[6] * axis.inner(n) > 0.0 ) axis.multiply_col(0,-1.0);
  axis.el(3) = 1.0;

  o.el(0) = a[11];
  o.el(1) = a[12];
  o.el(2) = a[13];
  o.el(3) = 1.0;

  axis = A * axis;
  n = A * n;
  o = A * o;

  a[5] = acos(axis.el(2));
  a[4] = acos(axis.el(0) / sin(a[5]));
  if (axis.el(0) * cos(a[4]) * sin(a[5]) < 0.0) a[5] = -a[5];
    else if (axis.el(1) * sin(a[4]) * sin(a[5]) < 0.0) a[4] = -a[4];

  a[3] = acos(n.el(2));
  a[2] = acos(n.el(0) / sin(a[3]));
  if (n.el(0) * cos(a[2]) * sin(a[3]) < 0.0) a[3] = -a[3];
    else if (n.el(1) * sin(a[2]) * sin(a[3]) < 0.0) a[2] = -a[2];
  
  a[11] = o.el(0);
  a[12] = o.el(1);
  a[13] = o.el(2);

  }

void cone::translate(vector3D& v)
{ a[11] += v.el(0);
  a[12] += v.el(1);
  a[13] += v.el(2);
  }

int cone::model_ok()
{ vector3D axis,n;
  double kot,rad,rad1; 

  axis.el(0) = cos(a[4]) * sin(a[5]);
  axis.el(1) = sin(a[4]) * sin(a[5]);
  axis.el(2) = cos(a[5]);
  if (a[6] > 0.0) axis.multiply_col(0,-1.0);
  n.el(0) = cos(a[2]) * sin(a[3]);
  n.el(1) = sin(a[2]) * sin(a[3]);
  n.el(2) = cos(a[3]);

  kot = acos(n.inner(axis));
  rad = fabs(zmin / tan(kot));
  rad1 = fabs(zmax / tan(kot));

  return (rad1 < max_radius && rad < max_radius);
  }


void cone::print()
{ vector3D axis,n;
  double kot,rad,rad1; 

  axis.el(0) = cos(a[4]) * sin(a[5]);
  axis.el(1) = sin(a[4]) * sin(a[5]);
  axis.el(2) = cos(a[5]);
  if (a[6] > 0.0) axis.multiply_col(0,-1.0);
  n.el(0) = cos(a[2]) * sin(a[3]);
  n.el(1) = sin(a[2]) * sin(a[3]);
  n.el(2) = cos(a[3]);

  kot = acos(n.inner(axis));
  rad = fabs(zmin / tan(kot));
  rad1 = fabs(zmax / tan(kot));

  printf("Cone: axis:(%lf,%lf,%lf)  n:(%lf,%lf,%lf)  radius1:%lf  radius2:%lf\n",axis.el(0),
   axis.el(1), axis.el(2), n.el(0), n.el(1), n.el(2), rad, rad1); 
 

}



void cone::get_normal_vector(double &x, double &y, double &z)
{ x = cos(a[4]) * sin(a[5]);
  y = sin(a[4]) * sin(a[5]);
  z = cos(a[5]);
  if (a[6] > 0.0)
  { x = -x;
    y = -y;
    z = -z;
    }
  }

double cone::max_radius = 1e20;
