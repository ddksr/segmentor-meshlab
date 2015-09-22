#include "matrix.h"

#include <math.h>
#include <stdlib.h>

extern "C" {
  #include "rcad_ndm_ut.h"
  void djacobi(double **a, int n, double *d, double **v, int *nrot);
  void deigsrt(double *d, double **v, int n);
}


// vo consists of original points
// vt consists of transformed points
// R will be a rotation matrix, so that || vt - R vo || will be smallest
// as possible 

void rotation(double *vo, double *vt, int n, symatrix &R)
{ int i,j,k;
  symatrix A(4),AT(4),C(4),W(3),W2(3),RS(3);
  double **B,**U,eval[5],gamma;
//  char s[100];    

  B = dmatrix(1,4,1,4);
  U = dmatrix(1,4,1,4);

  for (i = 1; i <= 4; i++)
    for (j = 1; j <= 4; j++) 
      B[i][j] = 0.0;
  
  for (i = 0; i < 4; i++) A.el(i,i) = 0.0;
   
  for (k = 0; k < n; k++)
  { 
    A.el(0,1) = vo[3*k+2] + vt[3*k+2];
    A.el(0,2) = -vo[3*k+1] - vt[3*k+1];
    A.el(1,2) = vo[3*k] + vt[3*k];
    for (i = 0; i < 3; i++) A.el(i,3) = vo[3*k+i] - vt[3*k+i];
    
    for (i = 1; i < 4; i++)
      for (j = 0; j < i; j++)
        A.el(i,j) = -A.el(j,i);

    AT = A.transpose();

    C = A * AT;

    for (i = 0; i < 4; i++)
      for (j = 0; j < 4; j++)
      { B[i+1][j+1] += C.el(i,j);
         }
    }

  djacobi(B,4,eval,U,&i);
  deigsrt(eval,U,4);

  R.identity();
  RS.identity();

  gamma = 2 * acos(U[4][4]);
  
  for (i = 0; i < 3; i++) W.el(i,i) = 0.0;
  W.el(0,1) = U[3][4] / sin(gamma/2);
  W.el(0,2) = -U[2][4] / sin(gamma/2);
  W.el(1,2) = U[1][4] / sin(gamma/2);
  for (i = 1; i < 3; i++)  
    for (j = 0; j < i; j++)
      W.el(i,j) = - W.el(j,i);
 
  W2 = W * W;

  for (i = 0; i < 3; i++)
  { W.multiply_col(i,sin(gamma));
    W2.multiply_col(i,1-cos(gamma));
    }  

  RS += (W + W2);

  RS = RS.transpose();

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      R.el(i,j) = RS.el(i,j);

  free_dmatrix(B,1,4,1,4);
  free_dmatrix(U,1,4,1,4);
  }

// vo consists of original points
// vt consists of transformed points
// t will be the vector, so that || vt - vo - t ||
// will be as smallest as possible

void translation(double *vo, double *vt, int n, vector3D &t)
{ int i;
  double invn = 1.0/n;
  t.el(0) = t.el(1) = t.el(2) = 0.0;

  for (i = 0; i < n; i++)
  { t.el(0) += (vt[3*i] - vo[3*i]);
    t.el(1) += (vt[3*i+1] - vo[3*i+1]);
    t.el(2) += (vt[3*i+2] - vo[3*i+2]);
    }

  t.multiply_col(0,invn);
  }

// vo consists of original points
// vt consists of transformed points
// A will be a rotation matrix, so that || vt - A vo || will be smallest
// as possible 

void transformation(double *vo, double *vt, int n, symatrix &A)
{ int i;
  double Mxx,Mxy,Mxz,Myy,Myz,Mzz,Mx,My,Mz;
  double Txx,Tyx,Tzx,Tx,Txy,Tyy,Tzy,Ty,Txz,Tyz,Tzz,Tz;
  double x,y,z,x1,y1,z1;
  symatrix M(4),MI(4);
  class matrix T(4,3),Ra(4,3),Tt(3,4);

  Mxx = Mxy = Mxz = Myy = Myz = Mzz = Mx = My = Mz = 0.0;
  Txx = Tyx = Tzx = Tx = Txy = Tyy = Tzy = Ty = Txz = Tyz = Tzz = Tz = 0.0;

  for (i = 0; i < n; i++)
  { x = vo[3*i];
    y = vo[3*i+1];
    z = vo[3*i+2];
    x1 = vt[3*i];
    y1 = vt[3*i+1];
    z1 = vt[3*i+2];

    Mxx += x*x;
    Mxy += x*y;
    Mxz += x*z;
    Myy += y*y;
    Myz += y*z;
    Mzz += z*z;
    Mx += x;
    My += y;
    Mz += z;
    Txx += x * x1;
    Txy += x * y1;
    Txz += x * z1;
    Tyx += y * x1;
    Tyy += y * y1;
    Tyz += y * z1;
    Tzx += z * x1;
    Tzy += z * y1;
    Tzz += z * z1;
    Tx += x1;
    Ty += y1;
    Tz += z1;
    }

  M.el(0,0) = Mxx;
  M.el(1,0) = M.el(0,1) = Mxy;
  M.el(2,0) = M.el(0,2) = Mxz;
  M.el(1,1) = M.el(1,1) = Myy;
  M.el(3,0) = M.el(0,3) = Mx;
  M.el(1,2) = M.el(2,1) = Myz;
  M.el(1,3) = M.el(3,1) = My;
  M.el(2,2) = Mzz;
  M.el(2,3) = M.el(3,2) = Mz;
  M.el(3,3) = n;
  
  T.el(0,0) = Txx;
  T.el(1,0) = Tyx;
  T.el(2,0) = Tzx;
  T.el(3,0) = Tx;
  T.el(0,1) = Txy;
  T.el(1,1) = Tyy;
  T.el(2,1) = Tzy;
  T.el(3,1) = Ty;
  T.el(0,2) = Txz;
  T.el(1,2) = Tyz;
  T.el(2,2) = Tzz;
  T.el(3,2) = Tz;  

  MI = M.inverse(i);

  Ra = MI * T;
  Tt = Ra.transpose();

  A.set_submatrix(Tt,0,0);
  A.el(3,0) = A.el(3,1) = A.el(3,2) = 0.0;
  A.el(3,3) = 1.0;

  }

void p_rotation(double *vo, double *vt, int n, symatrix &A)
{ double *von,*vtn,det[2];
  int i,j,k;
  symatrix R(3);
  vector3D vr;

  von = new double[3*n];
  vtn = new double[3*n];

  for (i = k = 0; i < n; i++)
  { 
    for (j = 0; j < 3; j++)
    { von[3*k+j] = vo[3*i+j];
      vtn[3*k+j] = vt[3*i+j];
      }
    det[0] = sqrt(von[3*k]*von[3*k]+von[3*k+1]*von[3*k+1]+von[3*k+2]*von[3*k+2]);
    det[1] = sqrt(vtn[3*k]*vtn[3*k]+vtn[3*k+1]*vtn[3*k+1]+vtn[3*k+2]*vtn[3*k+2]);
    if (det[0] > 1e-6 && det[1] > 1e-6)
   { for (j = 0; j < 3; j++)
      { von[3*k+j] /= det[0];
         vtn[3*k+j] /= det[1];
        }
      k++;
      }
   }

  rotation(von,vtn,k,R);
  A.identity();
  A.set_submatrix(R,0,0);

  delete [] von;
  delete [] vtn;

  }


/*
void p_rotation(double *vo, double *vt, int n, symatrix &A)
{ jpvector *r1, *r2;
   jpmatrix *Vr1, *Vr2, R(3,3), Rp(3,3), Rm(3,3);
   int i,j,it;
   double eps2;

   r1 = new jpvector[n];
   r2 = new jpvector[n];
   Vr1 = new jpmatrix[n];
   Vr2 = new jpmatrix[n];

   for (i = 0; i < n; i++)
   { r1[i].init(3);
      r2[i].init(3);
      Vr1[i].init(3,3);
      Vr2[i].init(3,3);
 
       r1[i][0] = vo[3*i];  r2[i][0] = vt[3*i];
       r1[i][1] = vo[3*i+1];  r2[i][1] = vt[3*i+1];
       r1[i][2] = vo[3*i+2];  r2[i][2] = vt[3*i+2];

       Vr1[i][0][0] = Vr1[i][1][1] = Vr1[i][2][2] = Vr2[i][0][0] = Vr2[i][1][1] = Vr2[i][2][2] = 1.0;
       Vr1[i][0][1] = Vr1[i][1][1] = Vr1[i][0][2] = Vr1[i][2][0] = Vr1[i][1][2] = Vr1[i][2][1] = 0.0;
       Vr2[i][0][1] = Vr2[i][1][1] = Vr2[i][0][2] = Vr2[i][2][0] = Vr2[i][1][2] = Vr2[i][2][1] = 0.0;
       }

   it = rotest(r2,Vr2,r1,Vr1,n,R,Rp,Rm,eps2, 0.05, 50);

   printf("\nIterations:%d  eps2:%lf", it,eps2);

   for (i = 0; i < 3; i++)
      for (j = 0; j < 3; j++) A.el(i,j) = R[j][i];

   delete [] r1;
   delete [] r2;
   delete [] Vr1;
   delete [] Vr2;

   }
*/
/*
void main()
{ double vt[30],vo[30];
  symatrix R(3),Ro(3);
  int i,j,k,n = 9;
   
  double d = 100.0;
  
    
  R.el(0,0) = R.el(1,1) = R.el(1,0) = 0.70711;
  R.el(0,1) = -0.70711;
  R.el(2,2) = 1;
  R.el(2,0) = R.el(2,1) = R.el(0,2) = R.el(1,2) = 0;

  vo[0][0] = 1;
  vo[0][1] = 0;
  vo[0][2] = 0;
  vo[1][0] = 0.70711;
  vo[1][1] = 0.70711;
  vo[1][2] = 0;
  vo[2][0] = 0;
  vo[2][1] = 1;
  vo[2][2] = 0;
  vo[3][0] = -0.70711;
  vo[3][1] = 0.70711;
  vo[3][2] = 0;

  vt[0][0] = 0.70711;
  vt[0][1] = 0.70711;
  vt[0][2] = 0;
  vt[1][0] = 0;
  vt[1][1] = 1;
  vt[1][2] = 0;
  vt[2][0] = -0.70711;
  vt[2][1] = 0.70711;
  vt[2][2] = 0;
  vt[3][0] = -1;
  vt[3][1] = 0;
  vt[3][2] = 0;

  
  for (i = 0; i < n; i++)  
  { for (j = 0; j < 2; j++)
    { vo[i][j] = 1.2*(rand()%11)/10.0 - 0.7;
      }

    vo[i][2] = sqrt(1 - vo[i][0] * vo[i][0] - vo[i][1] * vo[i][1]);

    
    vt[i][0] = R.el(0,0) * vo[i][0] + R.el(0,1) * vo[i][1] + R.el(0,2) * vo[i][2] + (rand()%10)/100.0;
    vt[i][1] = R.el(1,0) * vo[i][0] + R.el(1,1) * vo[i][1] + R.el(1,2) * vo[i][2] + (rand()%10)/100.0;
    vt[i][2] = R.el(2,0) * vo[i][0] + R.el(2,1) * vo[i][1] + R.el(2,2) * vo[i][2] + (rand()%10)/100.0;

    }
 
  vo[0] = -0.005860;
  vo[1] = -0.726328;
  vo[2] = 0.687323;
  vo[3] = -0.005581;
  vo[4] = 0.724788;
  vo[5] = 0.688949;
  vo[6] = -0.737292;
  vo[7] = -0.494566;
  vo[8] = 0.460222;

  vt[0] = 0.517604;
  vt[1] = -0.705697;
  vt[2] = 0.483816;
  vt[3] = 0.524241;
  vt[4] = 0.707226;
  vt[5] = 0.474345;
  vt[6] = -0.176405;
  vt[7] = -0.520416;
  vt[8] = 0.835493;

  rotation(vo,vt,2,Ro);

  cout << R << "\n" << Ro << "\n" << Ro.determinant();
 
  vector3D v,vout;

  v.el(0) = 1;
  v.el(1) = 2;
  v.el(2) = 3;

  vo[0] = 1.0;
  vo[1] = 0.0;
  vo[2] = -1.0;
  vo[3] = 2;
  vo[4] = 3;
  vo[5] = 0;

  for (i = 0; i < 6; i++) vt[i] = vo[i] + v.el(i % 3);
  
  translation(vo,vt,2,vout);

  cout << v << "\n" << vout;

  }

*/
