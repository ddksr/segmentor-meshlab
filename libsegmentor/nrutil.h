#ifndef LIBSEGMENTOR_NRUTIL
#define LIBSEGMENTOR_NRUTIL 1

void pmem();
float *vector(int nl, int nh);
float **matrix(int, int, int, int);
float **convert_matrix(float *, int, int, int, int);
double *dvector(int, int);
double **dmatrix(int, int, int, int);
int *ivector(int nl, int nh);
int **imatrix(int, int, int, int);
float **submatrix(float **, int, int, int, int, int, int);
void free_vector(float *, int, int);
void free_dvector(double *, int, int);
void free_ivector(int *, int, int);
void free_matrix(float **, int, int, int, int);
void free_dmatrix(double **, int, int, int, int);
void free_imatrix(int **, int, int, int, int);
void free_submatrix(float **, int, int, int, int);
void free_convert_matrix(float **, int, int, int, int);
void nrerror(char error_text[]);

#endif
