extern void estimate(struct vect *list, int no, double t[4][4]);
extern void convert(struct vect *list, int no, double T[4][4], double *a1, double *a2, double *a3, double *e1, double *e2, double *px, double *py, double *pz, double *fi, double *theta, double *psi );
extern int recover(struct vect *list, int no, double *a1, double *a2, double *a3, double *e1, double *e2, double *px, double *py, double *pz, double *fi, double *theta, double *psi );
extern int recover_search(struct vect *list, int no, double *a1, double *a2, double *a3, double *e1, double *e2, double *px, double *py, double *pz, double *fi, double *theta, double *psi, double *kx, double *ky, int rtype);

struct vect {
  double x;
  double y;
  double z;
};

#define RECOVER_SQ 1
#define RECOVER_SQ_TAPERING 2
