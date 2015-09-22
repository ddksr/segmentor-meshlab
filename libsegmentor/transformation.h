#ifndef LIBSEGMENTOR_TRANSFORMATION
#define LIBSEGMENTOR_TRANSFORMATION 1

#include "matrix.h"

void rotation(double *vo, double *vt, int n, symatrix& R);
void translation(double *vo, double *vt, int n, vector3D& t);
void transformation(double *vo, double *vt, int n, symatrix& A);
void p_rotation(double *vo, double *vt, int n, symatrix& A);

#endif
