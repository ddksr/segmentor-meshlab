/*----------------------------------------------------------------------------------*/
/*                                                                                  */
/* function to recover sq parameters of non-deformed sq from points                 */
/*                                                                                  */
/*----------------------------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define RAND_MAX        2147483647
#define mma 20

#include "deftype.h"
#define TRUE 1
#define FALSE 0

#define RECOVER_SQ 1
#define RECOVER_SQ_TAPERING 2
#define RECOVER_SQ_BENDING 3
#define RECOVER_SQ_GLOBAL 4
#define RECOVER_SQ_SYM_TAPERING 5
#define RECOVER_ASQ 6

/* #define PI 3.141529 */

struct vect {
  double x;
  double y;
  double z;
};



/* matherr handling */

/*
int matherr(exc)
struct exception *exc;
{
  fprintf(stderr, "%s(%f,%f)=%f\n", exc->name, exc->arg1, exc->arg2, exc->retval);
  return 0;

}

*/

/* -- begin of estimate.c -------------------------------------------------------------- */

estimate(list, no, t)
struct vect *list;
int no;
double t[4][4];
{
	double	a,      b,      c,
                x,	y,	z,
		xx,	yy,	zz,
		cxy,	cyz,	cxz,
		avgx,	avgy,	avgz,
                xmin,   ymin,   zmin,
                xmax,   ymax,   zmax,
		count, max_eigen;

	int	i,	j,	k;

	double	d[3],   aa[3][3],m_brain[3][3];
        double	m_inverse[4][4], determinant();

	int     nrot;

        x = y = z = 0; xx = yy = zz = cxy = cyz = cxz = 0; count = 0;
        xmin = ymin = zmin = 1000000.00;
        xmax = ymax = zmax = -1000000.00;

/* 	find	center of mass, cov, var and min and max extensions */

	for(i = 0; i < no; i++) {
	  a = list[i].x;
	  b = list[i].y;
	  c = list[i].z;
	  

         		x += a;
			y += b;
			z += c;

			xx += a*a;
			yy += b*b;
			zz += c*c;
			
			cxz += a*c;
			cyz += b*c;
			cxy += a*b;

			count++;
	}

/* 	center of mass coordinates, covariances and variances */

        if (count < 1.0){
          printf ("there are no range points in the standard input\n");
	}
	avgx = x/count;
	avgy = y/count;
	avgz = z/count;

	xx  /= count;
	yy  /= count;
	zz  /= count;

	cxz /= count;
	cyz /= count;
	cxy /= count;

        xx  -= (avgx * avgx);
        yy  -= (avgy * avgy);
        zz  -= (avgz * avgz);

        cxy -= (avgx * avgy);
        cxz -= (avgx * avgz);
        cyz -= (avgy * avgz);

/* 	arrange the numbers into matrix */

	aa[0][0] =  xx;   aa[0][1] = cxy; aa[0][2] = cxz;
	aa[1][0] = cxy;   aa[1][1] =  yy; aa[1][2] = cyz;
	aa[2][0] = cxz;   aa[2][1] = cyz; aa[2][2] =  zz;

/*
	printf("covariance matrix\n");
	print_matrix(aa);
*/

/* 	compute eigenvectors, eigenvalues */
	eigen (aa, 3, d, m_brain, &nrot);
/*
	printf("eigenvectors\n");
	print_matrix(m_brain);
	printf("eigenvalues (%15.8f,%15.8f,%15.8f)\n", d[0], d[1], d[2]);
	printf("\n");
        printf("stevilo rotacij: %d \n", nrot);
*/

/* along eigen vector with max eigen value is the smallest moment of inertia */

        max_eigen = -1000000000.0;

        if(d[0] > max_eigen) 
         { max_eigen = d[0]; i = 0;
	 }
        if(d[1] > max_eigen)
         { max_eigen = d[1]; i = 1;
	 }
        if(d[2] > max_eigen)
         { max_eigen = d[2]; i = 2;
	 }

        for (j = 0; j < 3; j++)
         { for (k = 0; k < 3; k++) m_inverse[j][k] = m_brain[j][k];
	 }
        for (j = 0; j < 3; j++)  
         { m_inverse[j][2] = m_brain[j][i];
           m_inverse[j][i] = m_brain[j][2];
	 }
        if (i != 2) for (j = 0; j < 3; j++) m_inverse[j][i] = -m_inverse[j][i];

        for (j = 0; j < 3; j++)
         { for (k = 0; k < 3; k++) m_brain[j][k] = m_inverse[j][k];
	 }

/* 	new axes must have roughly the same direction as old
   	for (i=0; i<3; i++){
	    if (m_brain[i][i] < 0){
                printf ("Eigenvector %d sign changed\n", i);
		for (j=0; j<3; j++){
		    m_brain[j][i] = - m_brain[j][i];
		}
	    }
	}
*/

	if( determinant(m_brain) < 0)
	    printf("determinant not positive\n");

/* 	print center of mass, object extensions, standard deviations */

/*	printf("center of mass: (%15.8f,%15.8f,%15.8f)\n", avgx, avgy, avgz);
*/

/* 	construct frame mat. in homogenuos coord. */

	for (i=0; i<3; i++)
	    for (j=0; j<3; j++){
		t[i][j] = m_brain[i][j];
	    }

	t[0][3] = avgx;
	t[1][3] = avgy;
	t[2][3] = avgz;
	t[3][3] =  1;

	/* Franc has forgotten the last row of the matrix ! A. Jaklic */

	t[3][0] = t[3][1] = t[3][2] = 0.0;

/* 	write out to standard output frame  matrix */

/*

	for (i=0; i<4; i++){
	   for (j=0; j<4; j++)
	   	printf("%18.10lf", t[i][j]);
	   printf("\n");
	}

*/

} /* end estimate */


/* Function: compute determinant of 3 x 3 matrix */

double determinant(matrix)
double	matrix[3][3];
{
double	d1, d2;

	d1 = matrix[0][0]*matrix[1][1]*matrix[2][2] +
	     matrix[0][1]*matrix[1][2]*matrix[2][0] +
	     matrix[0][2]*matrix[1][0]*matrix[2][1];
	d2 = matrix[0][2]*matrix[1][1]*matrix[2][0] +
	     matrix[0][0]*matrix[1][2]*matrix[2][1] +
	     matrix[0][1]*matrix[1][0]*matrix[2][2];

return(d1-d2);

} /* end determinant */


/* print matrix function */

print_matrix(matrix)
double	matrix[3][3];
{
int	i, j;

  for (i=0; i<3; i++){
    for (j=0; j<3; j++){
      printf("%15.8f", matrix[i][j]);
    }
    printf("\n");
  }
  printf("\n");

} /* end print_matrix */


/* Function for computing eigenvalues and eigenvectors, no sorting

            call: eigen (a, n, d, v, &nr);

            Input: a : pointer to the input array of type double
                   n : input (output) array size: n x n (int)

            Output:
                   d : pointer to the output array holding eigenvalues (double)
                   v : pointer to the output array holding eigenvectors (double)
                       v[i][j] column eigenvectors j, j = 0, 1, ..., n-1    
                  nr : number of iterations performed (int).

*/

int eigen( a, n,  d,  v, nr)
double *a;
int     n;
double *d;
double *v;
int    *nr;

{

int   nrot;

int   j, iq, ip, i;

double tresh,  theta,  tau,  t,  sm,  s,  h,  g,  c;

double b[256], z[256];


  for (ip = 0; ip < n; ip++){

    for (iq = 0; iq < n; iq++){

      v[ip * n + iq] = 0.0;
    }
    v[ip * n  + ip] = 1.0;
  }

  for (ip = 0; ip < n; ip++){

    b[ip] = a[ip * n + ip];
    d[ip] = b[ip];
    z[ip] = 0.0;
  }

  nrot = 0;

  for (i = 1; i < 50; i++){

    sm = 0.0;

    for ( ip = 0; ip < n-1; ip++){

      for ( iq = ip+1; iq <  n; iq++){

        sm += fabs( a[ip * n + iq]);

      }

    }

    if (sm == 0.0){
       *nr = nrot;
       return 0;
    }

    if (i < 6){

      tresh = (double)0.2 * sm / (double)(n * n);
    }

    else{

      tresh = 0.0;

    }

    for (ip = 0; ip < n-1; ip++){

      for (iq = ip+1; iq < n; iq++){

        g = 100.0 * fabs(a[ip * n + iq]);

        if ((i > 6) && ((fabs(d[ip]) + g) == fabs(d[ip]))
                    && ((fabs(d[iq]) + g) == fabs(d[iq]))){
           a[ip * n + iq] = 0.0;
        }

        else if (fabs(a[ip * n + iq]) > tresh){

          h = d[iq] - d[ip];

          if ((fabs(h) + g) == fabs(h)){

            t = a[ip * n + iq] / h;
	  }

          else{
            theta = 0.5 * h / a[ip * n + iq];
            t     = 1.0 /(fabs(theta) + sqrt(1.0 + theta * theta));
            if (theta < 0.0)  t = -t;
	  }

          c   = 1.0 / sqrt(1 + t * t);
          s   = t * c;
          tau = s / (1.0 + c);
          h   = t * a[ip * n + iq];

          z[ip] = z[ip] - h;
          z[iq] = z[iq] + h;
          d[ip] = d[ip] - h;
          d[iq] = d[iq] + h;

          a[ip * n + iq] = 0.0;

          for (j = 0; j <= ip-1; j++){

            g = a[j * n + ip];
            h = a[j * n + iq];
            a[j * n + ip] = g - s * (h + g * tau);
            a[j * n + iq] = h + s * (g - h * tau);
	  }

          for (j = ip+1; j <= iq-1; j++){

            g = a[ip * n + j];
            h = a[j  * n + iq];
            a[ip * n + j ] = g - s * (h + g * tau);
            a[j  * n + iq] = h + s * (g - h * tau);
	  }

          for (j = iq+1; j < n; j++){
            g = a[ip * n + j];
            h = a[iq * n + j];
            a[ip * n + j] = g - s * (h + g * tau);
            a[iq * n + j] = h + s * (g - h * tau);
	  }

          for (j = 0; j < n; j++){
            g = v[j * n + ip];
            h = v[j * n + iq];
            v[j * n + ip] = g - s * (h + g * tau);
            v[j * n + iq] = h + s * (g - h * tau);
	  }
          nrot++;

	}
      }
    }

    for (ip = 0; ip < n; ip++){

      b[ip] = b[ip] + z[ip];
      d[ip] = b[ip];
      z[ip] = 0.0;
    }
  }
  printf ("routine JACOBI: 50 iterations should not happen\n");
  *nr = nrot;
  return -1;

} /* end eigen */

/* -- end of estimate.c -------------------------------------------------------------- */


/* -- begin of convert.c ------------------------------------------------------------- */


/* given list of points of no elements and T compute initial parameters */

convert(list, no, T, a1, a2, a3, e1, e2, px, py, pz, fi, theta, psi )
struct vect list[];
int no;
double T[4][4];
double *a1,*a2, *a3, *e1, *e2, *px, *py, *pz, *fi, *theta, *psi;

{ int i, j, k;

  double r1, r2, r3, vector[4], result[4], a, b, c,
         xmin, xmax, ymin, ymax, zmin, zmax, TIN[4][4], ta1, ta2, ta3;

  char ans[10];

  double cos();
  double sin();
  double atan2();

/* 
reading input 
*/

/*
  now this is argumen T

  scanf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
                      &T[0][0], &T[0][1], &T[0][2], &T[0][3],
                      &T[1][0], &T[1][1], &T[1][2], &T[1][3],
                      &T[2][0], &T[2][1], &T[2][2], &T[2][3],
                      &T[3][0], &T[3][1], &T[3][2], &T[3][3]);
*/

  xmin = ymin = zmin = 100000.00;
  xmax = ymax = zmax = -100000.00;

  matrix_inverse(T, TIN);

  for(i = 0; i < no; i++) {
    vector[0] = list[i].x;
    vector[1] = list[i].y;
    vector[2] = list[i].z;
    vector[3] = 1;

     matrix_mult(vector, TIN, result);
     a = result[0]; b = result[1]; c = result[2];

     if(xmin > a)  xmin = a;
     
     if(xmax < a)  xmax = a;

     if(ymin > b)  ymin = b;
     
     if(ymax < b)  ymax = b;

     if(zmin > c)  zmin = c;
     
     if(zmax < c)  zmax = c;
   }


  r1 = atan2(T[1][2], T[0][2]);
  r1 = r1 + PI;

  r2 = atan2(cos(r1) * T[0][2] + sin(r1) * T[1][2], T[2][2]);

  r3 = atan2(-sin(r1) * T[0][0] + cos(r1) * T[1][0],
             -sin(r1) * T[0][1] + cos(r1) * T[1][1]);


  /* assign the return values */

  *fi    = r1;
  *theta = r2;
  *psi   = r3;

  /*  original code

  *px    = T[0][3];   
  *py    = T[1][3];
  *pz    = T[2][3];

  */

  /* to counter offset effect A. Jaklic 17.4.99 */ 


  *px    = T[0][3] + (xmax + xmin)/2;   
  *py    = T[1][3] + (ymax + ymin)/2;
  *pz    = T[2][3] + (zmax + zmin)/2;

 
  ta1    = (xmax - xmin)/2;
  ta2    = (ymax - ymin)/2;
  ta3    = (zmax - zmin)/2;


  ta1 = (ta1 > 1.0) ? ta1 : 1.0;        /* limit minimal SQ size A. Jaklic */
  ta2 = (ta2 > 1.0) ? ta2 : 1.0;
  ta3 = (ta3 > 1.0) ? ta3 : 1.0;

  *a1 = ta1;
  *a2 = ta2;
  *a3 = ta3;

  *e1    = 1;
  *e2    = 1;

/*

  printf("%c\n", 'n');

  printf("%f %f %f \n %f %f %f \n", r1, r2, r3, 
          T[0][3], T[1][3], T[2][3]);
  printf("%f %f %f\n", (xmax - xmin)/2, 
                       (ymax - ymin)/2, 
                       (zmax - zmin)/2);
  printf("\n\n");

  printf("%d %d \n", 1, 1);
  printf("%d %d \n", 0, 0);
  printf("%d %d %d %f %d \n", 0, -500, 500, 0.00001, 0);

  //printf("%d %d %f\n", 0, 0, 0.00001);

  printf("%d\n", 11); // 11 - undeformed sq, 15 - deformed sq 
  printf("%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", 
           0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 
              11, 12, 16, 17, 13, 14, 17, 18, 19, 20);

*/

}

/************************************************************************/
matrix_inverse(T, TIN)
double T[4][4], TIN[4][4];
{
 TIN[0][0] = T[0][0];
 TIN[0][1] = T[1][0];
 TIN[0][2] = T[2][0];
 TIN[0][3] = -(T[0][0]*T[0][3] + T[1][0]*T[1][3] + T[2][0]*T[2][3]);
 
 TIN[1][0] = T[0][1];
 TIN[1][1] = T[1][1];
 TIN[1][2] = T[2][1];
 TIN[1][3] = -(T[0][1]*T[0][3] + T[1][1]*T[1][3] + T[2][1]*T[2][3]);

 TIN[2][0] = T[0][2];
 TIN[2][1] = T[1][2];
 TIN[2][2] = T[2][2];
 TIN[2][3] = -(T[0][2]*T[0][3] + T[1][2]*T[1][3] + T[2][2]*T[2][3]);

 TIN[3][0] = 0;
 TIN[3][1] = 0;
 TIN[3][2] = 0;
 TIN[3][3] = 1;
}

/****************************************************************************
function for multiplying matrix with a vector
****************************************************************************/
 
matrix_mult(vector, matrix, result)
double vector[4], matrix[4][4], result[4];
{
  result[0] = matrix[0][0] * vector[0] +
              matrix[0][1] * vector[1] +
              matrix[0][2] * vector[2] +
              matrix[0][3] * vector[3];

  result[1] = matrix[1][0] * vector[0] +
              matrix[1][1] * vector[1] +
              matrix[1][2] * vector[2] +
              matrix[1][3] * vector[3];

  result[2] = matrix[2][0] * vector[0] +
              matrix[2][1] * vector[1] +
              matrix[2][2] * vector[2] +
              matrix[2][3] * vector[3];

  result[3] = matrix[3][0] * vector[0] +
              matrix[3][1] * vector[1] +
              matrix[3][2] * vector[2] +
              matrix[3][3] * vector[3];
}

/* -- begin of rec10.c -------------------------------------------------------------- */

/*---------------------------------------------------------------------------*/
/*									     */
/* rec10.c - least-square fitting					     */
/*									     */	
/*---------------------------------------------------------------------------*/

/* parameters convention

a[0] - a1
a[1] - a2
a[2] - a3
a[3] - e1
a[4] - e2
a[5] - fi
a[6] - theta
a[7] - psi
a[8] - px
a[9] - py
a[10]- pz
a[11]- kx
a[12]- ky
a[13]- center of bending   	;this parameter is allways fixed
a[14]- zmin			;this parameter is allways fixed
a[15]- zmax			;this parameter is allways fixed
a[16]- k
a[17]- alpha 

   end of convention */

int gliset;
double glochisq, gloldm;
glmma glatry, glbeta;
int i_am_in_trouble;

double mrqmin_init(), mrqmin(), mrqcof(), funcs(), sqr();


recover(list, no, a1, a2, a3, e1, e2, px, py, pz, fi, theta, psi)
struct vect list[];
int no;
double *a1, *a2, *a3, *e1, *e2, *px, *py, *pz, *fi, *theta, *psi;
{
 extern int gliset;
 extern double glochisq, gloldm;
 extern glmma glatry, glbeta;

 int npt, i, k, mfit, n_model_acc, iter;
 /* double mrqmin_init(), mrqmin(), mrqcof(), funcs(), sqr(); */
 double alamda, old_chisq;

 glmma a;
 glndata F, sig;
 glndata2 xw;
 gllista lista;
 glcovar covar, alpha;
 unsigned seed;

 FILE *fopen(), *fp;

 seed = 1955;
 srand(seed);

 gliset = 0;
 gloldm = -1.0;
 i_am_in_trouble = FALSE;


/*

 scanf("%*s");
 scanf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
         &a[5], &a[6], &a[7],
         &a[8], &a[9], &a[10],
         &a[0], &a[1], &a[2],
         &a[3], &a[4],
         &a[11], &a[12],
         &a[13], &a[14], &a[15], &a[16], &a[17]);
*/

 a[0]  = *a1;
 a[1]  = *a2;
 a[2]  = *a3;
 a[3]  = *e1;
 a[4]  = *e2;
 a[5]  = *fi;
 a[6]  = *theta;
 a[7]  = *psi;
 a[8]  = *px;
 a[9]  = *py;
 a[10] = *pz;

/* default values od deformation parameters to produce a nondeformed SQ */

 a[11] = 0;           /* kx */
 a[12] = 0;           /* ky */
 a[13] = 0;           /* center of bending */
 a[14] = -500;        /* z min */
 a[15] = 500;         /* z max */
 a[16] = 0.0000010;   /* k or 1/k ?? */
 a[17] = 0;           /* alpha angle */
 a[18] = 0;           /* tapering */
 a[19] = 1;           /* asq k */

 mfit = 11;  /* only first 11 parameters changed to minimize the function */

/*  scanf("%d", &mfit); */


 lista[0] = 0;
 lista[1] = 1;
 lista[2] = 2;
 lista[3] = 3;
 lista[4] = 4;
 lista[5] = 5;
 lista[6] = 6;
 lista[7] = 7;
 lista[8] = 8;
 lista[9] = 9;
 lista[10] = 10;
 lista[11] = 11;
 lista[12] = 12;
 lista[13] = 16;
 lista[14] = 17;
 lista[15] = 13;
 lista[16] = 14;
 lista[17] = 15;
 lista[18] = 18;
 lista[19] = 19;
 lista[20] = 20;



 /* scanf("%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",
         &lista[0], &lista[1], &lista[2], 
         &lista[3], &lista[4], 
         &lista[5], &lista[6], &lista[7], 
         &lista[8], &lista[9], &lista[10],
         &lista[11], &lista[12],
         &lista[13], &lista[14], &lista[15], &lista[16], &lista[17]);
 */


/**************************************************************************
read surface points:
*/
  i = 0;

  for (i = 0; i < no; i++) {
    xw[0][i] = list[i].x;
    xw[1][i] = list[i].y;
    xw[2][i] = list[i].z;
    F[i] = 0;
    sig[i] = 1;
  }

   

     npt = no;

     n_model_acc = npt;
     glochisq = 10e31;
     old_chisq = 1.0;
     iter = 0;

 /* printf("No. of points read [max. %d]: %d\n", arraysize, npt); */

     for (k = 0; k < 50; k++) {
       if (k == 0) 
          alamda = mrqmin_init(xw, F, sig, npt, a, lista, mfit, alpha, mma, &n_model_acc);
       else 
          alamda = mrqmin(xw, F, sig, npt, a,lista, mfit, covar, alpha, alamda, &n_model_acc);

       if (i_am_in_trouble)
	 return(FALSE);

        if(old_chisq  != glochisq) { 
	  old_chisq = glochisq;
	  
	  /*   print_to_file(iter, a, glochisq); */

          /*

           printf("Iteration # %d    alamda = %lf         chisq = %lf\n",
                    iter, alamda, glochisq);
           printf("a[0]        a[1]        a[2]        a[3]        a[4]\n");
           printf("%lf   %lf   %lf   %lf   %lf\n",
                  a[0], a[1], a[2], a[3], a[4]);
           printf("a[5]        a[6]        a[7]\n");
           printf("%lf   %lf   %lf\n", a[5], a[6], a[7]);
           printf("a[8]        a[9]        a[10]       a[11]       a[12]\n");
           printf("%lf   %lf   %lf   %lf   %lf\n",
                  a[8], a[9], a[10], a[11], a[12]);
           printf("a[13]       a[14]       a[15]       a[16]       a[17]\n");
           printf("%lf   %lf   %lf   %lf   %lf\n",
                  a[13], a[14], a[15], a[16], a[17]);

	   */

           iter = iter + 1;
	}
     }
	 
   
	 
	
  /* set return values */

   *a1 = a[0];
   *a2 = a[1];
   *a3 = a[2];
   *e1 = a[3];
   *e2 = a[4];
   *fi = a[5];
   *theta = a[6];
   *psi = a[7];
   *px = a[8];
   *py = a[9];
   *pz = a[10];

   return(TRUE);
}
	 
recover_search(list, no, a1, a2, a3, e1, e2, px, py, pz, fi, theta, psi, kx, ky, bk, ba, rtype)
struct vect list[];
int no;
double *a1, *a2, *a3, *e1, *e2, *px, *py, *pz, *fi, *theta, *psi, *kx, *ky, *bk, *ba;
int rtype;
{
 extern int gliset;
 extern double glochisq, gloldm;
 extern glmma glatry, glbeta;

 int npt, i, k, mfit, n_model_acc, iter;
 /* double mrqmin_init(), mrqmin(), mrqcof(), funcs(), sqr(); */
 double alamda, old_chisq;

 glmma a;
 glndata F, sig;
 glndata2 xw;
 gllista lista;
 glcovar covar, alpha;
 unsigned seed;

 FILE *fopen(), *fp;

 seed = 1955;
 srand(seed);

 gliset = 0;
 gloldm = -1.0;
 i_am_in_trouble = FALSE;


/*

 scanf("%*s");
 scanf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
         &a[5], &a[6], &a[7],
         &a[8], &a[9], &a[10],
         &a[0], &a[1], &a[2],
         &a[3], &a[4],
         &a[11], &a[12],
         &a[13], &a[14], &a[15], &a[16], &a[17]);
*/

 a[0]  = *a1;
 a[1]  = *a2;
 a[2]  = *a3;
 a[3]  = *e1;
 a[4]  = *e2;
 a[5]  = *fi;
 a[6]  = *theta;
 a[7]  = *psi;
 a[8]  = *px;
 a[9]  = *py;
 a[10] = *pz;

/* default values od deformation parameters to produce a nondeformed SQ */

 a[11] = *kx;           /* kx */
 a[12] = *ky;           /* ky */
 a[13] = 0;           /* center of bending */
 a[14] = -500;        /* z min */
 a[15] = 500;         /* z max */
 a[16] = *bk;   /* k or 1/k ?? */
 a[17] = *ba;           /* alpha angle */
 a[18] = 0;          /* sym tapering */
 a[19] = 0;          /* asq k */

 lista[0] = 5;
 lista[1] = 6;
 lista[2] = 7;
 lista[3] = 8;
 lista[4] = 9;
 lista[5] = 10;
 lista[6] = 11;
 lista[7] = 12;
 lista[8] = 16;
 lista[9] = 17;
 lista[10] = 18;
 lista[11] = 19;
 lista[12] = 0;
 lista[13] = 1;
 lista[14] = 2;
 lista[15] = 3;
 lista[16] = 4;
 lista[17] = 13;
 lista[18] = 14;
 lista[19] = 15;
 lista[20] = 20;

 switch(rtype) {
 case RECOVER_SQ_TAPERING:
   mfit = 8;
   break;
 case RECOVER_SQ_BENDING:
   mfit = 8;
   lista[6] = 16;
   lista[7] = 17;
   lista[8] = 11;
   lista[9] = 12;
   break;
 case RECOVER_SQ_SYM_TAPERING:
   lista[6] = 18;
   lista[7] = 19;
   lista[10] = 11;
   lista[11] = 19;
   mfit = 7;
   break;
 case RECOVER_ASQ:
   lista[6] = 18;
   lista[7] = 19;
   lista[10] = 11;
   lista[11] = 19;
   mfit = 8;
   break;
 case RECOVER_SQ:
 default:
   mfit = 6;
 }

 /* scanf("%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",
         &lista[0], &lista[1], &lista[2], 
         &lista[3], &lista[4], 
         &lista[5], &lista[6], &lista[7], 
         &lista[8], &lista[9], &lista[10],
         &lista[11], &lista[12],
         &lista[13], &lista[14], &lista[15], &lista[16], &lista[17]);
 */


/**************************************************************************
read surface points:
*/
  i = 0;

  for (i = 0; i < no; i++) {
    xw[0][i] = list[i].x;
    xw[1][i] = list[i].y;
    xw[2][i] = list[i].z;
    F[i] = 0;
    sig[i] = 1;
  }

   

     npt = no;

     n_model_acc = npt;
     glochisq = 10e31;
     old_chisq = 1.0;
     iter = 0;

 /* printf("No. of points read [max. %d]: %d\n", arraysize, npt); */

     for (k = 0; k < 50; k++) {
       if (k == 0) 
          alamda = mrqmin_init(xw, F, sig, npt, a, lista, mfit, alpha, mma, &n_model_acc);
       else 
          alamda = mrqmin(xw, F, sig, npt, a,lista, mfit, covar, alpha, alamda, &n_model_acc);

       if (i_am_in_trouble)
	 return(FALSE);

        if(old_chisq  != glochisq) { 
	  old_chisq = glochisq;
	  
	  /*   print_to_file(iter, a, glochisq); */

          /*

           printf("Iteration # %d    alamda = %lf         chisq = %lf\n",
                    iter, alamda, glochisq);
           printf("a[0]        a[1]        a[2]        a[3]        a[4]\n");
           printf("%lf   %lf   %lf   %lf   %lf\n",
                  a[0], a[1], a[2], a[3], a[4]);
           printf("a[5]        a[6]        a[7]\n");
           printf("%lf   %lf   %lf\n", a[5], a[6], a[7]);
           printf("a[8]        a[9]        a[10]       a[11]       a[12]\n");
           printf("%lf   %lf   %lf   %lf   %lf\n",
                  a[8], a[9], a[10], a[11], a[12]);
           printf("a[13]       a[14]       a[15]       a[16]       a[17]\n");
           printf("%lf   %lf   %lf   %lf   %lf\n",
                  a[13], a[14], a[15], a[16], a[17]);

	   */

           iter = iter + 1;
	}
     }
	 
   
	 
	
  /* set return values */

   *a1 = a[0];
   *a2 = a[1];
   *a3 = a[2];
   *e1 = a[3];
   *e2 = a[4];
   *fi = a[5];
   *theta = a[6];
   *psi = a[7];
   *px = a[8];
   *py = a[9];
   *pz = a[10];
   *kx = a[11];
   *ky = a[12];

   return(TRUE);
}	

recover2(list, no, a1, a2, a3, e1, e2, px, py, pz, fi, theta, psi, kx, ky, bk, ba, symt, asqk, rtype)
struct vect list[];
int no;
double *a1, *a2, *a3, *e1, *e2, *px, *py, *pz, *fi, *theta, *psi, *kx, *ky, *bk, *ba, *symt, *asqk;
int rtype;
{
 extern int gliset;
 extern double glochisq, gloldm;
 extern glmma glatry, glbeta;

 int npt, i, k, mfit, n_model_acc, iter;
 /* double mrqmin_init(), mrqmin(), mrqcof(), funcs(), sqr(); */
 double alamda, old_chisq;

 glmma a;
 glndata F, sig;
 glndata2 xw;
 gllista lista;
 glcovar covar, alpha;
 unsigned seed;

 FILE *fopen(), *fp;

 seed = 1955;
 srand(seed);

 gliset = 0;
 gloldm = -1.0;
 i_am_in_trouble = FALSE;


/*

 scanf("%*s");
 scanf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
         &a[5], &a[6], &a[7],
         &a[8], &a[9], &a[10],
         &a[0], &a[1], &a[2],
         &a[3], &a[4],
         &a[11], &a[12],
         &a[13], &a[14], &a[15], &a[16], &a[17]);
*/

 a[0]  = *a1;
 a[1]  = *a2;
 a[2]  = *a3;
 a[3]  = *e1;
 a[4]  = *e2;
 a[5]  = *fi;
 a[6]  = *theta;
 a[7]  = *psi;
 a[8]  = *px;
 a[9]  = *py;
 a[10] = *pz;

/* default values od deformation parameters to produce a nondeformed SQ */

 a[11] = *kx;           /* kx: 0 */
 a[12] = *ky;           /* ky: 0 */
 a[13] = 0;           /* center of bending */
 a[14] = -500;        /* z min */
 a[15] = 500;         /* z max */
 a[16] = *bk;   /* k or 1/k ?? : 0.0000010  */
 a[17] = *ba;           /* alpha angle : 0 */
 a[18] = *symt;   /* sym tapering */
 a[19] = *asqk;   /* asq k */

 mfit = 11;  /* only first 11 parameters changed to minimize the function */

/*  scanf("%d", &mfit); */


 lista[0] = 0;
 lista[1] = 1;
 lista[2] = 2;
 lista[3] = 3;
 lista[4] = 4;
 lista[5] = 5;
 lista[6] = 6;
 lista[7] = 7;
 lista[8] = 8;
 lista[9] = 9;
 lista[10] = 10;
 lista[11] = 11;
 lista[12] = 12;
 lista[13] = 16;
 lista[14] = 17;
 lista[15] = 18;
 lista[16] = 19;
 lista[17] = 15;
 lista[18] = 13;
 lista[19] = 14;
 lista[20] = 20;

 switch(rtype) {
 case RECOVER_SQ_TAPERING:
   mfit = 13;
   break;
 case RECOVER_SQ_BENDING:
   mfit = 13;
   lista[11] = 16;
   lista[12] = 17;
   lista[13] = 11;
   lista[14] = 12;
   break;
 case RECOVER_ASQ:
   mfit = 13;
   lista[11] = 18;
   lista[12] = 19;
   lista[15] = 11;
   lista[16] = 12; 
 case RECOVER_SQ_SYM_TAPERING:
   mfit = 12;
   lista[11] = 18;
   lista[12] = 19;
   lista[15] = 11;
   lista[16] = 12;   
   break;
 case RECOVER_SQ:
 default:
   mfit = 11;
 }


 /* scanf("%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",
         &lista[0], &lista[1], &lista[2], 
         &lista[3], &lista[4], 
         &lista[5], &lista[6], &lista[7], 
         &lista[8], &lista[9], &lista[10],
         &lista[11], &lista[12],
         &lista[13], &lista[14], &lista[15], &lista[16], &lista[17]);
 */


/**************************************************************************
read surface points:
*/
  i = 0;

  for (i = 0; i < no; i++) {
    xw[0][i] = list[i].x;
    xw[1][i] = list[i].y;
    xw[2][i] = list[i].z;
    F[i] = 0;
    sig[i] = 1;
  }

   

     npt = no;

     n_model_acc = npt;
     glochisq = 10e31;
     old_chisq = 1.0;
     iter = 0;

 /* printf("No. of points read [max. %d]: %d\n", arraysize, npt); */

     for (k = 0; k < 50; k++) {
       if (k == 0) 
          alamda = mrqmin_init(xw, F, sig, npt, a, lista, mfit, alpha, mma, &n_model_acc);
       else 
          alamda = mrqmin(xw, F, sig, npt, a,lista, mfit, covar, alpha, alamda, &n_model_acc);

       if (i_am_in_trouble)
	 return(FALSE);

        if(old_chisq  != glochisq) { 
	  old_chisq = glochisq;
	  
	  /*   print_to_file(iter, a, glochisq); */

          /*

           printf("Iteration # %d    alamda = %lf         chisq = %lf\n",
                    iter, alamda, glochisq);
           printf("a[0]        a[1]        a[2]        a[3]        a[4]\n");
           printf("%lf   %lf   %lf   %lf   %lf\n",
                  a[0], a[1], a[2], a[3], a[4]);
           printf("a[5]        a[6]        a[7]\n");
           printf("%lf   %lf   %lf\n", a[5], a[6], a[7]);
           printf("a[8]        a[9]        a[10]       a[11]       a[12]\n");
           printf("%lf   %lf   %lf   %lf   %lf\n",
                  a[8], a[9], a[10], a[11], a[12]);
           printf("a[13]       a[14]       a[15]       a[16]       a[17]\n");
           printf("%lf   %lf   %lf   %lf   %lf\n",
                  a[13], a[14], a[15], a[16], a[17]);

	   */

           iter = iter + 1;
	}
     }
	 
   
	 
	
  /* set return values */

   *a1 = a[0];
   *a2 = a[1];
   *a3 = a[2];
   *e1 = a[3];
   *e2 = a[4];
   *fi = a[5];
   *theta = a[6];
   *psi = a[7];
   *px = a[8];
   *py = a[9];
   *pz = a[10];
   *kx = a[11];
   *ky = a[12];
   *bk = a[16];
   *ba = a[17];
   *symt = a[18];
   *asqk = a[19];

   return(TRUE);
}

/***************************************************************************/
double mrqmin_init
  (x, F, sig, ndata, a, lista, mfit, alpha, nca, n_model_acc)
glndata2 x;
glndata F, sig;
int ndata;
glmma a;
gllista lista;
int mfit;
glcovar alpha;
int nca;
int *n_model_acc;
{ 
  int k, kk, j, ihit, n_model;
  double alamda, addnoise;

  kk = mfit;
  for(j = 0; j < nca; j++)
   { ihit = 0;
     for(k = 0; k < mfit; k++) 
      { if (lista[k] == j) ihit = ihit + 1;
      }

     if (ihit == 0)
      { lista[kk] = j; 
	kk = kk + 1;
      }
     else if (ihit > 1)
      { printf("pause 1 in routine MRQMIN\n");
        printf("Improper permutation in LISTA\n");
       
      }
   }

  if(kk != nca) 
   { printf("pause 2 in routine MRQMIN\n");
     printf("Improper permutation in LISTA\n");
   }

  alamda = 1;
  addnoise = 1e20;
  glochisq = mrqcof(x, F, sig, ndata, a, lista, mfit, alpha, glbeta, 
                    &n_model, n_model_acc, addnoise);
  *n_model_acc = n_model;
  return(alamda);
}

/***************************************************************************/
double mrqmin
   (x, F, sig, ndata, a, lista, mfit, covar, alpha, alamda, n_model_acc)
glndata2 x;
glndata F, sig;
int ndata;
glmma a;
gllista lista;
int mfit;
glcovar covar, alpha;
double alamda;
int *n_model_acc;
{ 
  double chisq;
  double poisson, addnoise;
  int rand();
  int k, j;
  glmma da;
  glcovar oneda;
  int n_model;

  for(j = 0; j < mfit; j++)
   { for(k = 0; k < mfit; k++) covar[j][k] = alpha[j][k];
     covar[j][j] = alpha[j][j] * (0.1 + alamda);
     oneda[j][0] = glbeta[j];
   }

  if (gaussj(covar, mfit, oneda) == FALSE) {
     i_am_in_trouble = TRUE;
     return(0);
  }
  for(j = 0; j < mfit; j++)
   { da[j] = oneda[j][0];
   }

  for (j = 0; j < mfit; j++)
   { if((lista[j] == 0) || (lista[j] == 1) || (lista[j] == 2))
      { if (a[lista[j]] + da[j] < 1.0) 
         { glatry[lista[j]] = 1.0;
	 /* printf("fixed a[0, 1, or 2]\n") */;
	 }
        else glatry[lista[j]] = a[lista[j]] + da[j];
      }
   /* omeji e1 in e2 na interval [0.1, 1.0] A. Jaklic julij 95, put back to 2.0 oct 96 to prevent minimization to get stuck*/
     else if((lista[j] == 3) || (lista[j] == 4))
      { if(a[lista[j]] + da[j] < 0.1) glatry[lista[j]] = 0.1;
        else if(a[lista[j]] + da[j] > 2.0) glatry[lista[j]] = 2.0;
        else glatry[lista[j]] = a[lista[j]] + da[j];
      }

     else if((lista[j] == 5) || (lista[j] == 6) || (lista[j] == 7))
      { glatry[lista[j]] = a[lista[j]] + da[j];
      }

     else if((lista[j] == 11) || (lista[j] == 12))
      { if(a[lista[j]] + da[j] > 1) glatry[lista[j]] = 1;
        else if(a[lista[j]] + da[j] < -1) glatry[lista[j]] = -1;
        else glatry[lista[j]] = a[lista[j]] + da[j];
      }

     else if(lista[j] == 13)
      { if(a[lista[j]] + da[j] < a[14] ||
           a[lista[j]] + da[j] > a[15]) glatry[lista[j]] = a[lista[j]];
        else glatry[lista[j]] = a[lista[j]] + da[j];
      }

     else if(lista[j] == 14)
      { if(a[lista[j]] + da[j] > a[13]) glatry[lista[j]] = a[lista[j]];
        else glatry[lista[j]] = a[lista[j]] + da[j];
      }

     else if(lista[j] == 15)
      { if(a[lista[j]] + da[j] < a[13]) glatry[lista[j]] = a[lista[j]];
        else glatry[lista[j]] = a[lista[j]] + da[j];
      }

     else if(lista[j] == 16)
      { if(fabs(a[lista[j]] + da[j]) > 1e9) glatry[lista[j]] = a[lista[j]];
        else if(a[lista[j]] + da[j] < 0.0) glatry[lista[j]] = a[lista[j]];
        else glatry[lista[j]] = a[lista[j]] + da[j];
      }

     else if(lista[j] == 17)
      {
        glatry[lista[j]] = a[lista[j]] + da[j];
        glatry[lista[j]] =
          glatry[lista[j]] - 2 * PI * (int)(glatry[lista[j]]/(2 * PI));
      }
	 else if(lista[j] == 18)
      { if(a[lista[j]] + da[j] > 1) glatry[lista[j]] = 1;
        else if(a[lista[j]] + da[j] < -1) glatry[lista[j]] = -1;
        else glatry[lista[j]] = a[lista[j]] + da[j];
      }
	 else if(lista[j] == 19) // amp
      {
        if(a[lista[j]] + da[j] < 0.0) glatry[lista[j]] = 0.0;
        else if(a[lista[j]] + da[j] > 1.0) glatry[lista[j]] = 1.0;
        else glatry[lista[j]] = a[lista[j]] + da[j];
      }
     else glatry[lista[j]] = a[lista[j]]+ da[j];
   }


 /*

  printf("------------------------------------------\n");
  printf(" Changes in parameters                    \n");
  printf("------------------------------------------\n");

  for(j = 0; j < mfit; j++)
   { printf("da[%2d] = %lf\n",lista[j], da[j]); 
   }

  */

  for (j = mfit; j < mma; j++) glatry[lista[j]] = a[lista[j]];

/***********************************************
ADDING NOISE
************************************************/

  poisson = ((double)rand())/RAND_MAX;  /* (4*32767.0); */
  addnoise = glochisq + (glochisq) * fabs(poisson)/2;

  chisq = mrqcof(x, F, sig, ndata, glatry, lista, mfit, covar, da,
                 &n_model, n_model_acc, addnoise);

  
  /*

  printf("\n");
  printf("           new_chisq: %lf           noise: %lf\n", chisq,poisson);
  printf("            glochisq: %lf\n", glochisq);
  printf("noise added glochisq: %lf\n", addnoise);

  */  

  if (chisq < addnoise)
   { if( chisq < addnoise) alamda = alamda/nu;
     glochisq = chisq;
     *n_model_acc = n_model;
     for(j = 0; j < mfit; j++)
      { glbeta[j] = da[j];
        a[lista[j]] = glatry[lista[j]];
        for(k = 0; k < mfit; k++) alpha[j][k] = covar[j][k];
      }
    /*  printf("\nACCEPT\n"); */
   }
  else
   { alamda = nu * alamda;
     chisq = glochisq;
     /* printf("\nREJECT\n"); */
   }

  return(alamda);
}

/*************************************************************************/
double mrqcof
    (x, F, sig, ndata, a, lista, mfit, alpha, beta, n_model, n_model_acc, 
     addnoise)
glndata2 x;
glndata F, sig;
int ndata;
glmma a;
gllista lista;
int mfit;
glcovar alpha;
glmma beta;
int *n_model;
int *n_model_acc;
double addnoise;
{
  int k, j, i;
  double Fmod, wt, sig2i, dF, chisq, treshold;
  glmma dFda;
  double rotx[11], roty[11], rotz[11];
  
  for(j= 0; j < mfit; j++)
   { for(k = 0; k <= j; k++) alpha[j][k] = 0.0;
     beta[j] = 0.0;
   }

  sig2i = 1.0/(sig[1] * sig[1]);
  precomp_rot(a, rotx, roty, rotz);

/*****************  SEGMENTATION ********************************/
/*
  treshold = 2 * sqrt(glochisq/sig2i);
*/
  treshold = 10e30;

  do
   { 
  /* printf("\n"); */

/*     printf("treshold = %lf     ", treshold);
*/
     chisq = 0.0;
     *n_model = 0;
     for(i = 0; i < ndata; i++)
      { Fmod = funcs(x[0][i], x[1][i], x[2][i], rotx, roty, rotz, a, dFda);
        dF = F[i] - Fmod;

        if (Fmod < treshold)
         { (*n_model) +=  1;
           for(j = 0; j < mfit; j++)
            { wt = dFda[lista[j]]*sig2i;
  
              for(k = 0; k <= j; k++)
                   alpha[j][k] = alpha[j][k] + wt * dFda[lista[k]];
              beta[j] = beta[j] + dF*wt;
   	    }
           chisq = chisq + dF*dF*sig2i;
         }
        else 
         {/* printf("point rejected: %f\n", Fmod);*/
         }
/*****  The iteration cannot be accepted  *********/

        if(chisq / (*n_model_acc) > addnoise) break;
      }
     treshold = treshold * 2;
/*
  printf("accepted points: %d    rejected points: %d \n",
                                          (*n_model), ndata - (*n_model));
*/
     if(chisq / (*n_model_acc) > addnoise) break;
   } 
  while ((double)(*n_model)/(*n_model_acc) < 0.95);

  for(j = 0; j < mfit; j++)
   { for(k = 0; k <= j ; k++) 
      { alpha[j][k] = alpha[j][k] * chisq;
      }
     beta[j] = beta[j] * chisq;
   }
 
  chisq = chisq/(*n_model);
/*
  printf("accepted points: %d    rejected points: %d \n",
                                          (*n_model), ndata - (*n_model));
*/
  for(j = 1; j < mfit; j++)
   { for(k = 0; k <= j-1; k++) 
      { alpha[k][j] = alpha[j][k];
      }
   }
  return(chisq);
}

/****************************************************************************/
int gaussj(a, n, b)
glcovar a, b;
int n;

{ double big, dum, pivinv;
  int i, icol, irow, j, k, l, ll, m;
  glnp indxc, indxr, ipiv;

  m = 1;

  for(j = 0; j < n; j++)
   { ipiv[j] = 0;
   }

  for(i = 0; i < n; i++)
   { big = 0.0;
     for(j = 0; j < n; j++)
      { if (ipiv[j] != 1)
         { for(k = 0; k < n ; k++)
            { if(ipiv[k] == 0)
               { if (fabs(a[j][k]) >= big)
                  { big = fabs(a[j][k]);
  		    irow = j;
	 	    icol = k;
		  }
	       }
	    else if (ipiv[k] > 1) {
                         printf("pause 1 in GAUSSJ - singular matrix\n");
			 return(FALSE);
			 /* exit(1); */
	                 }
	    }
	 }
      }

     ipiv[icol] = ipiv[icol] +1;
     if (irow != icol) 
      {	for( l = 0; l < n; l++)
         { dum = a[irow][l];
	   a[irow][l] = a[icol][l];
	   a[icol][l] = dum;
	 }
        for(l = 0; l < m; l++)
         { dum = b[irow][l];
	   b[irow][l] = b[icol][l];
	   b[icol][l] = dum;
	 }
      }

     indxr[i] = irow;
     indxc[i] = icol;
     if (a[icol][icol] == 0.0) printf("pause 2 in GAUSSJ - singular matrix\n");

     pivinv = 1.0/a[icol][icol];
     a[icol][icol] = 1.0;
     for(l = 0; l < n; l++) a[icol][l] = a[icol][l] *pivinv;
     for(l = 0; l < m; l ++) b[icol][l] = b[icol][l] * pivinv;
     for(ll = 0; ll < n; ll++)
      { if (ll != icol)
         { dum = a[ll][icol];
	   a[ll][icol] = 0.0;
	   for(l = 0; l < n; l++) a[ll][l] = a[ll][l] - a[icol][l]* dum;
	   for(l = 0; l < m; l++) b[ll][l] = b[ll][l] - b[icol][l]*dum;
	 }
      }
   }

  for(l = n-1; l >= 0; l--)
   { if (indxr[l] != indxc[l])
      { for( k = 0; k < n; k++)
         { dum = a[k][indxr[l]];
	   a[k][indxr[l]] = a[k][indxc[l]];
	   a[k][indxc[l]] = dum;
	 }
      }
   }
 return(TRUE);
}

/**************************************************************************/
precomp_rot(a, rotx, roty, rotz)
glnparam a;
double rotx[11], roty[11], rotz[11];
{ 
 rotx[0] = (cos(a[5])*cos(a[6])*cos(a[7]) - sin(a[5])*sin(a[7]));
 rotx[1] = (-cos(a[5])*cos(a[6])*sin(a[7]) - sin(a[5])*cos(a[7]));
 rotx[2] = (cos(a[5])*sin(a[6]));

 rotx[3] = (-sin(a[5])*cos(a[6])*cos(a[7]) - cos(a[5])*sin(a[7]));
 rotx[4] = (sin(a[5])*cos(a[6])*sin(a[7]) - cos(a[5])*cos(a[7]));
 rotx[5] = (-sin(a[5])*sin(a[6]));

 rotx[6] = (-cos(a[5])*sin(a[6])*cos(a[7]));
 rotx[7] = (cos(a[5])*sin(a[6])*sin(a[7]));
 rotx[8] = (cos(a[5])*cos(a[6]));

 rotx[9] = (-cos(a[5])*cos(a[6])*sin(a[7]) - sin(a[5])*cos(a[7]));
 rotx[10]= (-cos(a[5])*cos(a[6])*cos(a[7]) + sin(a[5])*sin(a[7]));

 roty[0] = (sin(a[5])*cos(a[6])*cos(a[7]) + cos(a[5])*sin(a[7]));
 roty[1] = (-sin(a[5])*cos(a[6])*sin(a[7]) + cos(a[5])*cos(a[7]));
 roty[2] = (sin(a[5])*sin(a[6]));

 roty[3] = (cos(a[5])*cos(a[6])*cos(a[7]) - sin(a[5])*sin(a[7]));
 roty[4] = (-cos(a[5])*cos(a[6])*sin(a[7]) - sin(a[5])*cos(a[7]));
 roty[5] = (cos(a[5])*sin(a[6]));
 
 roty[6] = (-sin(a[5])*sin(a[6])*cos(a[7]));
 roty[7] = (sin(a[5])*sin(a[6])*sin(a[7]));
 roty[8] = (sin(a[5])*cos(a[6]));

 roty[9] = (-sin(a[5])*cos(a[6])*sin(a[7]) + cos(a[5])*cos(a[7]));
 roty[10]= (-sin(a[5])*cos(a[6])*cos(a[7]) - cos(a[5])*sin(a[7]));

 rotz[0] = (-sin(a[6])*cos(a[7]));
 rotz[1] = (sin(a[6])*sin(a[7]));
 rotz[2] = (cos(a[6]));

 rotz[3] = 0;
 rotz[4] = 0;
 rotz[5] = 0;

 rotz[6] = (-cos(a[6])*cos(a[7]));
 rotz[7] = (cos(a[6])*sin(a[7]));
 rotz[8] = (-sin(a[6]));

 rotz[9] = (sin(a[6])*sin(a[7]));
 rotz[10] = (sin(a[6])*cos(a[7]));

}

double mylog(x)
double x;
{
  if (x > 0)
    return(log(x));
  else {
    printf("Log (%lf) = \n", x);
       return(0);
  }
    
}

double mypow(x, y)
double x, y;
{
  if (x >= 0)
    return(pow(x,y));
  else {
    printf("Pow (%lf, %lf) = %lf\n", x, y, pow(x,y));
 
    return(pow(fabs(x),y));
  }
    
}

/****************************************************************************/
double funcs(x, y, z, rotx, roty, rotz, a, dFda)
double x, y, z;
glnparam a, dFda;
double rotx[11], roty[11], rotz[11];
{ int i;
  double coefx, coefy, coefz, A, B, D, G, F, DA, DB, GD, Gz, R,
         M, N, nx, ny, nz, nz_p, N1, N2, SAB, CAB;
  glnparam dx, dy, dz, dxb, dyb, dzb, dA, dB, dD, dG, dR;
  double xb, yb, zb, koren;
  double O, Omin, Omax, zhat;

  // ASQ:
  glnparam dxa, dya, dza;
  double xa, ya, za;

  double mypow(), mylog();

 coefx = (rotx[0]*x+
          roty[0]*y + 
          rotz[0]*z - (a[8]*rotx[0] + a[9]*roty[0] + a[10]*rotz[0])
         );

 for(i = 0; i < 5; i++) dx[i] = 0;

 dx[5] = (rotx[3] * x +
          roty[3] * y +
          rotz[3] * z - (a[8]*rotx[3] + a[9]*roty[3] + a[10]*rotz[3])
         );

 dx[6] = (rotx[6] * x +
          roty[6] * y +
          rotz[6] * z - (a[8]*rotx[6] + a[9]*roty[6] + a[10]*rotz[6])
         );

 dx[7] = (rotx[9] * x +
          roty[9] * y +
          rotz[9] * z - (a[8]*rotx[9] + a[9]*roty[9] + a[10]*rotz[9])
         );

 dx[8] = -rotx[0];
 dx[9] = -roty[0];
 dx[10] = -rotz[0];

 for(i = 11; i < mma; i++) dx[i] = 0;

 coefy = (rotx[1]*x+
          roty[1]*y + 
          rotz[1]*z - (a[8]*rotx[1] + a[9]*roty[1] + a[10]*rotz[1])
         );

 for(i = 0; i < 5; i++) dy[i] = 0;

 dy[5] = (rotx[4] * x +
          roty[4] * y +
          rotz[4] * z - (a[8]*rotx[4] + a[9]*roty[4] + a[10]*rotz[4])
         );

 dy[6] = (rotx[7] * x +
          roty[7] * y +
          rotz[7] * z - (a[8]*rotx[7] + a[9]*roty[7] + a[10]*rotz[7])
         );

 dy[7] = (rotx[10] * x +
          roty[10] * y +
          rotz[10] * z - (a[8]*rotx[10] + a[9]*roty[10] + a[10]*rotz[10])
         );

 dy[8] = -rotx[1];
 dy[9] = -roty[1];
 dy[10] = -rotz[1];

 for(i = 11; i < mma; i++) dy[i] = 0;

 coefz = (rotx[2]*x+
          roty[2]*y + 
          rotz[2]*z - (a[8]*rotx[2] + a[9]*roty[2] + a[10]*rotz[2])
         );

 for(i = 0; i < 5; i++) dz[i] = 0;

 dz[5] = (rotx[5] * x +
          roty[5] * y +
          rotz[5] * z - (a[8]*rotx[5] + a[9]*roty[5] + a[10]*rotz[5])
         );

 dz[6] = (rotx[8] * x +
          roty[8] * y +
          rotz[8] * z - (a[8]*rotx[8] + a[9]*roty[8] + a[10]*rotz[8])
         );

 dz[7] = 0;

 dz[8] = -rotx[2];
 dz[9] = -roty[2];
 dz[10] = -rotz[2];

 for(i = 11; i < mma; i++) dz[i] = 0;

/********************************************/
/*
  xb = coefx;
  for (i = 0; i < mma; i++) dxb[i] = dx[i];

  yb = coefy;
  for (i = 0; i < mma; i++) dyb[i] = dy[i];

  zb = coefz;
  for (i = 0; i < mma; i++) dzb[i] = dz[i];
*/
/********************************************/

  koren = sqrt(sqr(coefx) + sqr(coefy));
  SAB = sin(a[17] - atan2(coefy, coefx));
  CAB = cos(a[17] - atan2(coefy, coefx));

  R = 1/(koren * CAB);

  for(i = 0; i < 5; i++) dR[i] = 0;
  for(i = 5; i < 11; i++)
   dR[i] = (-(coefx * dx[i] + coefy * dy[i]) / CAB -
            (SAB/sqr(CAB)) * (coefx * dy[i] - coefy * dx[i]))/
             (koren * sqr(koren));
  for (i = 11; i < 16; i++) dR[i] = 0;

  dR[17] = SAB/(koren * sqr(CAB));
  dR[18] = 0;

/**********************************************/

  koren = sqrt(sqr(coefz) + sqr(1/a[16] - 1/R));

  xb = coefx - (1/R - 1/a[16] + koren) * cos(a[17]);

  for(i = 0; i < 5; i++) dxb[i] = 0;
  for(i = 0; i < 11; i++)
   dxb[i] = dx[i] -
            (-dR[i]/sqr(R) + (1.0/koren) * 
                     (coefz * dz[i] + (1/a[16] - 1/R) * (dR[i]/sqr(R)))) * 
              cos(a[17]);
  for (i = 11; i < 16; i++) dxb[i] = 0;

  dxb[16] = cos(a[17]) * (-1/sqr(a[16]) + (1/a[16] - 1/R)/
                                          (koren * sqr(a[16])));

  dxb[17] = sin(a[17]) * (1/R - 1/a[16] + koren) +
            cos(a[17]) * dR[17] * (1/sqr(R) - (1/a[16] - 1/R)/(koren*sqr(R)));

  dxb[18] = 0;

/************************************************/

  yb = coefy - (1/R - 1/a[16] + koren) * sin(a[17]);

  for(i = 0; i < 5; i++) dyb[i] = 0;
  for(i = 0; i < 11; i++)
   dyb[i] = dy[i] -
            (-dR[i]/sqr(R) + (1.0/koren) * 
            (coefz * dz[i] + (1/a[16] - 1/R) * dR[i]/sqr(R))) * sin(a[17]);
  for (i = 11; i < 16; i++) dyb[i] = 0;

  dyb[16] = sin(a[17]) * (-1/sqr(a[16]) + (1/a[16] - 1/R)/
                                          (koren * sqr(a[16])));

  dyb[17] = -cos(a[17]) * (1/R - 1/a[16] + koren) -
            sin(a[17]) * dR[17] * (-1/sqr(R) + (1/a[16] - 1/R)/
                                               (koren* sqr(R)));
  dyb[18] = 0;

/***************************************************************/

  O = atan2(coefz, 1/a[16] - 1/R);

  zb = O/a[16];

  for (i = 0; i < 5; i++) dzb[i] = 0;
  for (i = 5; i < 11; i++) 
   dzb[i] = (1/(a[16]*(sqr(1/a[16] - 1/R) + sqr(coefz)))) *
            ((1/a[16] - 1/R) * dz[i] - coefz * dR[i]/sqr(R));
  for (i = 11; i < 16; i++) dzb[i] = 0;

  dzb[16] = -O/sqr(a[16]) +
            (1/(a[16] * sqr(a[16]))) * coefz/(sqr(coefz) + sqr(1/a[16] - 1/R));

  dzb[17] = 0;
  dzb[18] = 0;

/******************************************************************/

  // tmp pass
  for (i = 0; i < mma; i++) {
	dxa[i] = dxb[i];
	dya[i] = dyb[i];
	dza[i] = dzb[i];
  }
  xa = xb; ya = yb; za = zb;

  


/******************************************************************/  

/* OLD SIMPLE TAPERING 
 A = xa/(a[0]*(a[11]*za/a[2] + 1));

 dA[0] = -A/a[0];
 dA[1] = 0;
 dA[2] = (xa/a[0]) * (1.0/sqr(a[11]*za/a[2] + 1))*a[11]*za/sqr(a[2]);

 for(i = 3; i < 11; i++)
   dA[i] = dxa[i]/(a[0]*(a[11]*za/a[2] + 1)) - 
            xa * a[11] * dza[i]/(a[0]*sqr(a[11] * za/a[2] + 1)*a[2]);

   dA[11] = -xa * za/(a[0]*sqr(a[11]*za/a[2] + 1)*a[2]);

   dA[12] = 0;

 for(i = 13; i < mma; i++)
   dA[i] = dxa[i]/(a[0]*(a[11]*za/a[2] + 1)) - 
            xa * a[11] * dza[i]/(a[0]*sqr(a[11] * za/a[2] + 1)*a[2]);


 B = ya/(a[1]*(a[12]*za/a[2] + 1));

 dB[0] = 0;
 dB[1] = -B/a[1];
 dB[2] = (ya/a[1]) * (1.0/sqr(a[12]*za/a[2] + 1))*a[12]*za/sqr(a[2]);

 for(i = 3; i < 11; i++)
   dB[i] = dya[i]/(a[1]*(a[12]*za/a[2] + 1)) - 
            ya * a[12] * dza[i]/(a[1]*sqr(a[12] * za/a[2] + 1)*a[2]);

   dB[11] = 0;

   dB[12] = -ya * za/(a[1]*sqr(a[12]*za/a[2] + 1)*a[2]);

 for(i = 13; i < mma; i++)
   dB[i] = dya[i]/(a[1]*(a[12]*za/a[2] + 1)) - 
            ya * a[12] * dza[i]/(a[1]*sqr(a[12] * za/a[2] + 1)*a[2]);

*/

  /*********************************************/

  // Symetric tapering
  /*********************************************/

  
/*********************************************/
  // ASQ tapering

  A = xa/(a[0]*pow(pow(cos((a[19]*PI*za)/(2.*a[2])),2),0.5)*(1 + (a[11]*za)/a[2])*(1 + (a[18]*za)/a[2]));

  dA[0] = -(xa/(pow(a[0],2)*pow(pow(cos((a[19]*PI*za)/(2.*a[2])),2),0.5)*(1 + (a[11]*za)/a[2])*(1 + (a[18]*za)/a[2])));

  dA[1] = 0;

  dA[2] = (a[18]*xa*za)/(a[0]*pow(a[2],2)*pow(pow(cos((a[19]*PI*za)/(2.*a[2])),2),0.5)*(1 + (a[11]*za)/a[2])*pow(1 + (a[18]*za)/a[2],2)) + (a[11]*xa*za)/(a[0]*pow(a[2],2)*pow(pow(cos((a[19]*PI*za)/(2.*a[2])),2),0.5)*pow(1 + (a[11]*za)/a[2],2)*(1 + (a[18]*za)/a[2])) - (1.5707963267948966*a[19]*cos((a[19]*PI*za)/(2.*a[2]))*sin((a[19]*PI*za)/(2.*a[2]))*xa*za)/(a[0]*pow(a[2],2)*pow(pow(cos((a[19]*PI*za)/(2.*a[2])),2),1.5)*(1 + (a[11]*za)/a[2])*(1 + (a[18]*za)/a[2]));

  for (i = 3; i < 17; i++) {
	if (i == 11 || i == 12) { continue; }
	dA[i] = dxa[i]/(a[0]*pow(pow(cos((a[19]*PI*za)/(2.*a[2])),2),0.5)*(1 + (a[11]*za)/a[2])*(1 + (a[18]*za)/a[2])) - (a[18]*xa*dza[i])/(a[0]*a[2]*pow(pow(cos((a[19]*PI*za)/(2.*a[2])),2),0.5)*(1 + (a[11]*za)/a[2])*pow(1 + (a[18]*za)/a[2],2)) - (a[11]*xa*dza[i])/(a[0]*a[2]*pow(pow(cos((a[19]*PI*za)/(2.*a[2])),2),0.5)*pow(1 + (a[11]*za)/a[2],2)*(1 + (a[18]*za)/a[2])) + (1.5707963267948966*a[19]*cos((a[19]*PI*za)/(2.*a[2]))*sin((a[19]*PI*za)/(2.*a[2]))*xa*dza[i])/(a[0]*a[2]*pow(pow(cos((a[19]*PI*za)/(2.*a[2])),2),1.5)*(1 + (a[11]*za)/a[2])*(1 + (a[18]*za)/a[2]));
  }

  dA[11] = -((xa*za)/(a[0]*a[2]*pow(pow(cos((a[19]*PI*za)/(2.*a[2])),2),0.5)*pow(1 + (a[11]*za)/a[2],2)*(1 + (a[18]*za)/a[2])));
  dA[12] = 0;

  dA[18] = -((xa*za)/(a[0]*a[2]*pow(pow(cos((a[19]*PI*za)/(2.*a[2])),2),0.5)*(1 + (a[11]*za)/a[2])*pow(1 + (a[18]*za)/a[2],2)));

  dA[19] = (1.5707963267948966*cos((a[19]*PI*za)/(2.*a[2]))*sin((a[19]*PI*za)/(2.*a[2]))*xa*za)/(a[0]*a[2]*pow(pow(cos((a[19]*PI*za)/(2.*a[2])),2),1.5)*(1 + (a[11]*za)/a[2])*(1 + (a[18]*za)/a[2]));



  B = ya/(a[1]*pow(pow(cos((a[19]*PI*za)/(2.*a[2])),2),0.5)*(1 + (a[12]*za)/a[2])*(1 + (a[18]*za)/a[2]));
  dB[0] = 0;
  dB[1] = -(ya/(pow(a[1],2)*pow(pow(cos((a[19]*PI*za)/(2.*a[2])),2),0.5)*(1 + (a[12]*za)/a[2])*(1 + (a[18]*za)/a[2])));
  dB[2] = (a[18]*ya*za)/(a[1]*pow(a[2],2)*pow(pow(cos((a[19]*PI*za)/(2.*a[2])),2),0.5)*(1 + (a[12]*za)/a[2])*pow(1 + (a[18]*za)/a[2],2)) + (a[12]*ya*za)/(a[1]*pow(a[2],2)*pow(pow(cos((a[19]*PI*za)/(2.*a[2])),2),0.5)*pow(1 + (a[12]*za)/a[2],2)*(1 + (a[18]*za)/a[2])) - (1.5707963267948966*a[19]*cos((a[19]*PI*za)/(2.*a[2]))*sin((a[19]*PI*za)/(2.*a[2]))*ya*za)/(a[1]*pow(a[2],2)*pow(pow(cos((a[19]*PI*za)/(2.*a[2])),2),1.5)*(1 + (a[12]*za)/a[2])*(1 + (a[18]*za)/a[2]));
  for (i = 3; i < 17; i++) {
	if (i == 11 || i == 12) { continue; }
	dB[i] = dya[i]/(a[1]*pow(pow(cos((a[19]*PI*za)/(2.*a[2])),2),0.5)*(1 + (a[12]*za)/a[2])*(1 + (a[18]*za)/a[2])) - (a[18]*ya*dza[i])/(a[1]*a[2]*pow(pow(cos((a[19]*PI*za)/(2.*a[2])),2),0.5)*(1 + (a[12]*za)/a[2])*pow(1 + (a[18]*za)/a[2],2)) - (a[12]*ya*dza[i])/(a[1]*a[2]*pow(pow(cos((a[19]*PI*za)/(2.*a[2])),2),0.5)*pow(1 + (a[12]*za)/a[2],2)*(1 + (a[18]*za)/a[2])) + (1.5707963267948966*a[19]*cos((a[19]*PI*za)/(2.*a[2]))*sin((a[19]*PI*za)/(2.*a[2]))*ya*dza[i])/(a[1]*a[2]*pow(pow(cos((a[19]*PI*za)/(2.*a[2])),2),1.5)*(1 + (a[12]*za)/a[2])*(1 + (a[18]*za)/a[2]));
  }

  dB[11] = 0;
  dB[12] = -((ya*za)/(a[1]*a[2]*pow(pow(cos((a[19]*PI*za)/(2.*a[2])),2),0.5)*pow(1 + (a[12]*za)/a[2],2)*(1 + (a[18]*za)/a[2])));

  dB[18] = -((ya*za)/(a[1]*a[2]*pow(pow(cos((a[19]*PI*za)/(2.*a[2])),2),0.5)*(1 + (a[12]*za)/a[2])*pow(1 + (a[18]*za)/a[2],2)));
  dB[19] = (1.5707963267948966*cos((a[19]*PI*za)/(2.*a[2]))*sin((a[19]*PI*za)/(2.*a[2]))*ya*za)/(a[1]*a[2]*pow(pow(cos((a[19]*PI*za)/(2.*a[2])),2),1.5)*(1 + (a[12]*za)/a[2])*(1 + (a[18]*za)/a[2]));
  

/*********************************************/
 DA = mypow(sqr(A), 1.0/a[4]);   DB = mypow(sqr(B), 1.0/a[4]);
 D = DA + DB;
 
 for(i = 0; i < 3; i++)
   dD[i] = (2/a[4]) * (DA/A) * dA[i] +
            (2/a[4]) * (DB/B) * dB[i];

   dD[3] = 0;
   dD[4] = -DA * mylog(sqr(A))/sqr(a[4]) -
             DB * mylog(sqr(B))/sqr(a[4]);

 for(i = 5; i < mma; i++)
   dD[i] = (2/a[4]) * (DA/A) * dA[i] +
            (2/a[4]) * (DB/B) * dB[i];

/**********************************************/
 GD = mypow(D, a[4]/a[3]);   Gz = mypow(sqr(za/a[2]), 1.0/a[3]);
 G = GD + Gz;
 
 for(i = 0; i < 2; i++)
     dG[i] = (a[4]/a[3]) * (GD/D) * dD[i];

 dG[2] = (a[4]/a[3]) * (GD/D) * dD[2] -
         (2/a[3]) * (a[2]/za)* Gz * (za/sqr(a[2]));

 dG[3] = GD *
            (-a[4] * (1/sqr(a[3])) * mylog(D) + (a[4]/a[3])*(1/D) * dD[3]) -
            (1.0/sqr(a[3])) * Gz * mylog(sqr(za/a[2]));

 dG[4] = GD *
            ((1/a[3]) * mylog(D) + (a[4]/a[3])*(1/D) * dD[4]);

 for(i = 5; i < mma; i++)
     dG[i] = (a[4]/a[3]) * (GD/D) * dD[i] + 
             (2/a[3]) * (Gz/za) * dza[i];

/***********************************************/
 F = mypow(G, a[3]);
 for(i= 0; i < 3; i++)
   dFda[i] = a[3] * (F/G) * dG[i];

   dFda[3] = F * (mylog(G) + (a[3]/G) * dG[3]);

 for(i = 4; i < mma; i++)
   dFda[i] = a[3] * (F/G) * dG[i];

/****** correction of criteria function */
/*
        for(i = 0; i < mma; i++)
         dFda[i] = 0.5 * dFda[i] /sqrt(F);
        F = sqrt(fabs(F)) * F/fabs(F);
*/
        F = F - 1;

        dFda[0] = sqrt(a[1] * a[2]) * 
                  (F/(2*sqrt(a[0])) + sqrt(a[0]) * dFda[0]);
        dFda[1] = sqrt(a[0] * a[2]) * 
                  (F/(2*sqrt(a[1])) + sqrt(a[1]) * dFda[1]);
        dFda[2] = sqrt(a[0] * a[1]) * 
                  (F/(2*sqrt(a[2])) + sqrt(a[2]) * dFda[2]);

        for(i = 3; i < mma; i++) 
           dFda[i] = dFda[i] * sqrt(a[0] * a[1] * a[2]);

        F = sqrt(a[0] * a[1] * a[2]) * F;

/********** computing the normal in this point **************

        M = mypow(sqr(coefz), 1.0/a[3]);
        N1 = mypow(sqr(coefx), 1.0/a[4]);
        N2 = mypow(sqr(coefy), 1.0/a[4]);

        if(fabs(N2) < 10e-10) N2 = 10e-10;
        if(fabs(N1) < 10e-10) N1 = 10e-10;
        N = N1/N2;

        nx = (a[0]/coefx) * (N/(1 + N)) * (1 - M);
        ny = (a[1]/coefy) * (1.0/ (1 + N)) * (1 - M);
        nz = (a[2]/coefz) * M;

        nz_p = -sin(a[6])*cos(a[7]) * nx +
                sin(a[6])*sin(a[7]) * ny +
                cos(a[6])           * nz;

        if(nz_p < 0) 
         { printf("%5.0lf ", F);
           F = 0.1 * F;
	 }
*/
        if(F < 0) return(1 * F);
        else return(F);
}

/*************************************************************************/
double sqr(x)
double x;
{ return (x*x);
}

/*************************************************************************/
print_to_file(i, a, chisq)
int i;
glmma a;
double chisq;
{ char f[10];
  int one;
  double zero;
  FILE *fopen(), *fp;

  sprintf(f, "fit_%d", i);

  fp = fopen(f, "w");

/*  fprintf(fp, "n\n"); */
  fprintf(fp, " %lf %lf %lf \n", a[5], a[6], a[7]);
  fprintf(fp, " %lf %lf %lf \n", a[8], a[9], a[10]);
  fprintf(fp, " %lf %lf %lf \n", a[0], a[1], a[2]);
  zero = 0;
  fprintf(fp, " %lf %lf \n", zero, zero);
  fprintf(fp, " %lf %lf \n", a[3], a[4]);
  one = 1;
  fprintf(fp, " %d %d %d \n", 2, 4, 1); /* 2 4 1 */
  fprintf(fp, " %lf %lf\n", a[11], a[12]);
  fprintf(fp, " %lf %lf %lf %lf %lf\n", 0.0, -500.0, 500.0, 1.0/a[16], a[17]);
  fprintf(fp, "n\n");
  fprintf(fp, "%lf\n", chisq);

  fclose(fp);

  /* Create a file in format compatible with REC10.C input */

  sprintf(f, "iter_%d", i);

  fp = fopen(f, "w");

  fprintf(fp, "n\n");

  /* a[5] = fi, a[6] = theta, a[7] = psi */

  fprintf(fp, " %lf %lf %lf \n", a[5], a[6], a[7]);

  /* a[8] = px, a[9] = py, a[10] = pz */

  fprintf(fp, " %lf %lf %lf \n", a[8], a[9], a[10]);

  /* a[0] = a1, a[1] = a1, a[2] = a3 */

  fprintf(fp, " %lf %lf %lf \n", a[0], a[1], a[2]);

  /* a[3] = epsiilon1, a[4] = epsilon2 */

  fprintf(fp, " %lf %lf \n", a[3], a[4]);

  /* a[11] = Kx, a[12] = Ky */

  fprintf(fp, " %lf %lf\n", a[11], a[12]);

  /* a[13] = center of bending (constant !), a[14], a[15] = start and end
     of bend (zmin), (zmax) (constant !), a[16] = k, a[17] = alpha,
     radius of bending = 1/a[16] */

  fprintf(fp, " %lf %lf %lf %lf %lf\n", 0.0, -500.0, 500.0, a[16], a[17]);

  fclose(fp);

}




