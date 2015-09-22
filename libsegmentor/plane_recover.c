#include "rcad_ndm_ut.h"

#include <math.h>
#include <stdio.h>

#define SQR( x ) ((x)*(x))



#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);

void djacobi( a, n, d, v, nrot)
double  **a, *d, **v;
int     n, *nrot; 
{
	int j,iq,ip,i;
	double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;
	

	b=dvector(1,n);
	z=dvector(1,n);
	for (ip=1;ip<=n;ip++) {
		for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=1;ip<=n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	*nrot=0;
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0) {
			free_dvector(z,1,n);
			free_dvector(b,1,n);
			return;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++) {
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && fabs(d[ip])+g == fabs(d[ip])
					&& fabs(d[iq])+g == fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if (fabs(h)+g == fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=1;j<=ip-1;j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<=iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<=n;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=1;j<=n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
				}
			}
		}
		for (ip=1;ip<=n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
        for (i = 1; i <= n; i++)
	  { printf("\n");
            for (j = 1; j <= n; j++) printf("%lf ",a[i][j]);
	  }
	nrerror("Too many iterations in routine JACOBI");
      }

#undef ROTATE


void deigsrt( d, v, n)
double  *d, **v; 
int     n;
{
	int k,j,i;
	double p;

	for (i=1;i<n;i++) {
		p=d[k=i];
		for (j=i+1;j<=n;j++)
			if (d[j] >= p) p=d[k=j];
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			for (j=1;j<=n;j++) {
				p=v[j][i];
				v[j][i]=v[j][k];
				v[j][k]=p;
			}
		}
	}
}

int eigenfit( x,y, ndata, a, chisq)
int ndata;
double **x,y[],a[],*chisq;
{
  int       i;
  double    sum_x, sum_y, sum_z, norm;
  double    cen_x,cen_y,cen_z;
  double    sum_xx, sum_xy, sum_xz, sum_yy, sum_yz, sum_zz;
  double    xp,yp,zp;
  double    **mcm;                       /* matrix of central moments */
  double    *eigenvalue;                 /* eigenvalues of mcm */
  double    **eigenvector;          /* normalized eigenvectors of mcm */
  int       nrot; 


  mcm         = (double**) dmatrix(1,3,1,3);
  eigenvector = (double**) dmatrix(1,3,1,3);
  eigenvalue  = (double*)  dvector(1,3);

  sum_x = sum_y = sum_z = 0.0;
  sum_xx =sum_xy =sum_xz =sum_yy =sum_yz =sum_zz= 0.0;



  for(i=1; i<=ndata; i++) {
    xp = x[i][1];
    yp = x[i][2];
    zp = y[i];
    sum_x  += xp;
    sum_y  += yp;
    sum_z  += zp;
    sum_xx += xp*xp;
    sum_xy += xp*yp;
    sum_xz += xp*zp;
    sum_yy += yp*yp;
    sum_yz += yp*zp;
    sum_zz += zp*zp;
  }
    
  cen_x = sum_x/ndata;
  cen_y = sum_y/ndata;
  cen_z = sum_z/ndata;
 
  
  /*** fill in the mcm matrix ***/
  mcm[1][1] =              sum_xx - sum_x*sum_x/ndata;
  mcm[2][2] =              sum_yy - sum_y*sum_y/ndata;         
  mcm[3][3] =              sum_zz - sum_z*sum_z/ndata;
  mcm[1][2] = mcm[2][1] =  sum_xy - sum_x*sum_y/ndata;
  mcm[1][3] = mcm[3][1] =  sum_xz - sum_x*sum_z/ndata;
  mcm[2][3] = mcm[3][2] =  sum_yz - sum_y*sum_z/ndata;


  /*** jacobi eigenvector fit ***/
  djacobi(mcm,3,eigenvalue,eigenvector,&nrot);
  deigsrt(eigenvalue,eigenvector,3);            /* sort eigvec in desc. order */

  /*printf("\n\tEIGENVECTORS");
  for(i=1;i<=3;i++) { 
    printf("\n%lf  %lf %lf -> %lf",
	 eigenvector[1][i],eigenvector[2][i],eigenvector[3][i],eigenvalue[i]);
  }*/

  norm = sqrt(SQR(eigenvector[1][3])+SQR(eigenvector[2][3])+SQR(eigenvector[3][3]));
  eigenvector[1][3] /= norm;
  eigenvector[2][3] /= norm;
  eigenvector[3][3] /= norm;
 

  a[1] = eigenvector[1][3];
  a[2] = eigenvector[2][3];
  a[3] = eigenvector[3][3];

  /* ensure that the normal vector has positive z direction */
  if( a[3] < 0.0 )
    {
      a[1] = -a[1];
      a[2] = -a[2];
      a[3] = -a[3];
    }

  a[4] = -(cen_x*a[1] + cen_y*a[2] + cen_z*a[3]); 

  /* goodness of fit */

  *chisq = eigenvalue[3];



  free_dmatrix(mcm,1,3,1,3);
  free_dmatrix(eigenvector,1,3,1,3);
  free_dvector(eigenvalue,1,3);
  return(1);
}





