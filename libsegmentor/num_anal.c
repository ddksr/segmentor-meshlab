/* num_anal.c   c file containing NUmerical analysis (Numerical  Recipes)

Non-lin least squares fitting and associated  routines   */


#include <math.h>
#include "rcad_ndm_ut.h"

#define TINY 1.0e-20

#define SQR( x ) ((x)*(x))

void o_mrqmin(x,y,sig,ndata,a,ma,lista,mfit,covar,alpha,chisq,alamda, funcs)
double **x,y[],sig[],a[],**covar,**alpha,*chisq,*alamda;
int ndata,ma,lista[],mfit;
void (*funcs)(); 
{
	int k,kk,j,ihit;
	static double *da,*atry,**oneda,*beta,ochisq;
	
	void o_mrqcof(),o_gaussj(),o_covsrt(),nrerror(),free_dmatrix(),free_dvector();
    
	if (*alamda < 0.0) {
		oneda=dmatrix(1,mfit,1,1);
		atry=dvector(1,ma);
		da=dvector(1,ma);
		beta=dvector(1,ma);
		kk=mfit+1;
		for (j=1;j<=ma;j++) {
			ihit=0;
			for (k=1;k<=mfit;k++)
				if (lista[k] == j) ihit++;
			if (ihit == 0)
				lista[kk++]=j;
			else if (ihit > 1) nrerror("Bad LISTA permutation in MRQMIN-1");
		}
		if (kk != ma+1) nrerror("Bad LISTA permutation in MRQMIN-2");
		*alamda=0.001;
		/* mrqcof(x,y,sig,ndata,a,ma,lista,mfit,alpha,beta,chisq,funcs); */
		o_mrqcof(x,y,sig,ndata,a,ma,lista,mfit,alpha,beta,chisq, funcs); 
		
		ochisq=(*chisq);
	}
	for (j=1;j<=mfit;j++) {
		for (k=1;k<=mfit;k++) covar[j][k]=alpha[j][k];
		covar[j][j]=alpha[j][j]*(1.0+(*alamda));
		oneda[j][1]=beta[j];
	}
	o_gaussj(covar,mfit,oneda,1);
	for (j=1;j<=mfit;j++)
		da[j]=oneda[j][1];
	if (*alamda == 0.0) {
		o_covsrt(covar,ma,lista,mfit);
		free_dvector(beta,1,ma);
		free_dvector(da,1,ma);
		free_dvector(atry,1,ma);
		free_dmatrix(oneda,1,mfit,1,1);
		return;
	}
	for (j=1;j<=ma;j++) atry[j]=a[j];
	for (j=1;j<=mfit;j++)
		atry[lista[j]] = a[lista[j]]+da[j];
	/* mrqcof(x,y,sig,ndata,atry,ma,lista,mfit,covar,da,chisq,funcs); */
	o_mrqcof(x,y,sig,ndata,atry,ma,lista,mfit,covar,da,chisq, funcs); 
	if (*chisq < ochisq) {
		*alamda *= 0.1;
		ochisq=(*chisq);
		for (j=1;j<=mfit;j++) {
			for (k=1;k<=mfit;k++) alpha[j][k]=covar[j][k];
			beta[j]=da[j];
			a[lista[j]]=atry[lista[j]];
		}
	} else {
		*alamda *= 10.0;
		*chisq=ochisq;
	}
	return;
}


void o_mrqcof(x,y,sig,ndata,a,ma,lista,mfit,alpha,beta,chisq, funcs)
double **x,y[],sig[],a[],**alpha,beta[],*chisq;
int ndata,ma,lista[],mfit;
 void (*funcs)(); 	/* ANSI: void (*funcs)(float,float *,float *,float *,int); */
{
	int k,j,i;
	double ymod,wt,sig2i,dy,*dfda;
	
    
	dfda = (double *)dvector(1,ma);
	for (j=1;j<=mfit;j++) {
		for (k=1;k<=j;k++) alpha[j][k]=0.0;
		beta[j]=0.0;
	}
	*chisq=0.0;
	for (i=1;i<=ndata;i++) {
		(* funcs)(x[i][1],x[i][2],x[i][3],a,&ymod,dfda,ma);
		sig2i=1.0/(sig[i]*sig[i]);
		dy=y[i]-ymod;
		for (j=1;j<=mfit;j++) {
			wt=dfda[lista[j]]*sig2i;
			for (k=1;k<=j;k++)
				alpha[j][k] += wt*dfda[lista[k]];
			beta[j] += dy*wt;
		}
		(*chisq) += dy*dy*sig2i;
	}
	for (j=2;j<=mfit;j++)
		for (k=1;k<=j-1;k++) alpha[k][j]=alpha[j][k];
	free_dvector(dfda,1,ma);
	
	(*chisq) = (*chisq)/ (double) ndata;
}


    
#define SWAP(a,b) {float temp=(a);(a)=(b);(b)=temp;}
    
void o_gaussj(a,n,b,m)
double **a,**b;
int n,m;
{       int bad = 0;
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll,*ivector();
	double big,dum,pivinv;
	void nrerror(),free_ivector();
 
 icol = irow = 0;
 
	indxc=ivector(1,n);
	indxr=ivector(1,n);
	ipiv=ivector(1,n);
	for (j=1;j<=n;j++) ipiv[j]=0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if (ipiv[j] != 1)
				for (k=1;k<=n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} else if (ipiv[k] > 1) nrerror("GAUSSJ: Singular Matrix-1");
				}
  if (icol < 1 || icol > n || irow < 1 || irow > n)
		{ printf("\nBad Matrix %d %d %d",icol,irow,n);
    nrerror("GAUSSJ: Bad Matrix ");
                  bad = 1;
		  } else
		{ ++(ipiv[icol]);
		  if (irow != icol) {
			for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
			for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
		  }
		  indxr[i]=irow;
		  indxc[i]=icol;
		  if (a[icol][icol] == 0.0) 
                  { nrerror("GAUSSJ: Singular Matrix-2");
                    bad = 2;
		    } else
		  { pivinv=1.0/a[icol][icol];
		    a[icol][icol]=1.0;
		    for (l=1;l<=n;l++) a[icol][l] *= pivinv;
		    for (l=1;l<=m;l++) b[icol][l] *= pivinv;
		    for (ll=1;ll<=n;ll++)
			  if (ll != icol) {
				  dum=a[ll][icol];
				  a[ll][icol]=0.0;
				  for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
				  for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
			  }
		  }
		}	  
        }

        if (!bad)
	  for (l=n;l>=1;l--) {
		  if (indxr[l] != indxc[l])
			  for (k=1;k<=n;k++)
				  SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
}
    
#undef SWAP


void o_covsrt(covar,ma,lista,mfit)
double **covar;
int ma,lista[],mfit;
{
	int i,j;
	double swap;
    
	for (j=1;j<ma;j++)
		for (i=j+1;i<=ma;i++) covar[i][j]=0.0;
	for (i=1;i<mfit;i++)
		for (j=i+1;j<=mfit;j++) {
			if (lista[j] > lista[i])
				covar[lista[j]][lista[i]]=covar[i][j];
			else
				covar[lista[i]][lista[j]]=covar[i][j];
		}
	swap=covar[1][1];
	for (j=1;j<=ma;j++) {
		covar[1][j]=covar[j][j];
		covar[j][j]=0.0;
	}
	covar[lista[1]][lista[1]]=swap;
	for (j=2;j<=mfit;j++) covar[lista[j]][lista[j]]=covar[1][j];
	for (j=2;j<=ma;j++)
		for (i=1;i<=j-1;i++) covar[i][j]=covar[j][i];
}













