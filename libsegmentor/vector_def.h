/*  C++ interface to C vector manipulation routines */
#ifndef LIBSEGMENTOR_VECTORDEF
#define LIBSEGMENTOR_VECTORDEF


double dot ( double *, double *);
double modsq(double *);
void normalise(double *, double *);
void cross_prod( double *, double *, double *);
double dist(double , double , double , double , double , double );

void vector_mult_sc(double *, double);

#endif
