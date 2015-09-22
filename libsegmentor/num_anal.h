/*  Header file for non linear least sq estimations from Numerical Recipes */

#ifndef LIBSEGMENTOR_NUMERICAL_ANALYSIS
#define LIBSEGMENTOR_NUMERICAL_ANALYSIS 1

void o_mrqmin(double **x, double *y,double *sig,int ndata,double *a,int ma,int *lista,int mfit,double **covar,double **alpha,double *chisq,double *alamda, int (*funcs)(double,double ,double, double *,double *,double *,int));

void o_mrqcof(double **x, double *y,double *sig,int ndata,double **a,int ma,int *lista,int mfit,double **alpha, double *beta,double **chisq, int (*funcs)(double,double ,double, double *,double *,double *,int));



void o_gaussj(double **a,int n,double **b,int m);




void o_covsrt(double **covar,int ma,int lista,int mfit);

#endif

