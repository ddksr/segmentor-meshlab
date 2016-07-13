extern void estimate(struct vect *list, int no, double t[4][4]);
extern void convert(struct vect *list, int no, double T[4][4], double *a1, double *a2, double *a3, double *e1, double *e2, double *px, double *py, double *pz, double *fi, double *theta, double *psi );
extern int recover(struct vect *list, int no, double *a1, double *a2, double *a3, double *e1, double *e2, double *px, double *py, double *pz, double *fi, double *theta, double *psi );
extern int recover_search(struct vect *list, int no, double *a1, double *a2, double *a3, double *e1, double *e2, double *px, double *py, double *pz, double *fi, double *theta, double *psi, double *kx, double *ky, double *bk, double *ba, double *symt, double *asqk, int rtype);
extern int recover2(struct vect *list, int no, double *a1, double *a2, double *a3, double *e1, double *e2, double *px, double *py, double *pz, double *fi, double *theta, double *psi, double *kx, double *ky, double *bk, double *ba, double *symt, double *asqk, int rtype);

struct vect {
  double x;
  double y;
  double z;
};

#define RECOVER_SQ 1
#define RECOVER_SQ_TAPERING 2
#define RECOVER_SQ_BENDING 3
#define RECOVER_SQ_GLOBAL 4
#define RECOVER_SQ_SYM_TAPERING 5
#define RECOVER_ASQ 6
