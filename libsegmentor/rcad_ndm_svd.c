/* NDM_SVD.C */

/* 27-Aug-1991 LG */

/* Modifications:
		18-SEP-1991  AL
                12-NOV-1991  LG
                04-DEC-1991  LG
                12-MAR-1992  HT-KGY  -- see MINDOUBLE
                21-DEC-1993  LG-KZ Test if Nan in svdvmp
		21-APR-1997  LG ANSI/ not ANSI
*/


#include <stdio.h>
#include <math.h>
/* #include "g_memory.h" */
#include "rcad_ndm_ut.h"

#include "ndm_svd.h"
#include <float.h>

#ifndef HUGE
  #define HUGE  FLT_MAX
#endif

#ifdef __mips
#define  MINDOUBLE (1.0e-307)
#else
#define  MINDOUBLE (1.0/HUGE)
#endif

static int ITERNUM=100;

/* Singular value decomposition and least squares solution */
/* Press-Flannery-Teukolsky-Vetterling :Num. Rec. pp 65-71. */


typedef struct SVDDESC {
	double	*rv1;
	int	ncom;
} svddesc;


static svddesc glob = { 0, 0};

#ifdef ANSI
void init_svdcmp( int n, int flag, char ** regime)
#else
void init_svdcmp( n, flag, regime)
int	n, flag;
char	**regime;
#endif
/********************************/
/* n= number of columns!!! */
{
	svddesc * reg;

	reg = (svddesc * ) * regime;
	if (flag != 0) {
		if (reg == 0) {
			(*regime) = (char *) GIVMEM(doub_mh, sizeof(svddesc));
			reg = (svddesc * ) * regime;
			reg->rv1 = dvector(1, n);
			reg->ncom = n;
		} else if (n > reg->ncom)
			nrerror("Incompatible dimensions in svdcmp");
		glob = *reg;
	} else if ((flag == 0) && (reg != 0)) {
		free_dvector( reg->rv1, 1, n);
		RELMEM(doub_mh, (char *) reg);
		*regime = 0;
	}
}

#ifndef ANSI
void svbksb( u, w, v, m, n, b, x, rbuf)
/*************************************/

double	**u, w[], **v, b[], x[];
int	m, n;
char	**rbuf;
#else
void svbksb(double ** u, double w[], double **v, int m, int n,
            double b[], double x[], char **rbuf)
/*************************************/
#endif

{
	int	jj, j, i;
	double	s;

	init_svdcmp(n, 1, rbuf);
	for (j = 1; j <= n; j++) {
		s = 0.0;
		if (w[j]) {
			for (i = 1; i <= m; i++) 
				s += u[i][j] * b[i];
			s /= w[j];
		}
		(glob.rv1)[j] = s;
	}
	for (j = 1; j <= n; j++) {
		s = 0.0;
		for (jj = 1; jj <= n; jj++) 
			s += v[j][jj] * (glob.rv1)[jj];
		x[j] = s;
	}
}



static double	at, bt, ct;

#define PYTHAG(a,b)  ((at=fabs(a)) > (bt=fabs(b)) ? \
(ct=bt/at,at*sqrt(1.0+ct*ct)) : \
(bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)) : 0.0))


static double	maxarg1, maxarg2;

#define MAX(a,b)  (maxarg1=(a), maxarg2=(b), (maxarg1) > (maxarg2) ?\
	(maxarg1) : (maxarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


#ifndef ANSI
int svdcmp( a, m, n, w, v, rbuf)
/*******************************/

double	**a, *w, **v;
int	m, n;
char	**rbuf;
#else
int svdcmp( double **a, int m, int n, double *w,
             double **v, char **rbuf)
/*******************************/
#endif
{
	int	flag, i, its, j, jj, k, l, nm;
	double	c, f, h, s, x, y, z;
	double	anorm = 0.0, g = 0.0, scale = 0.0;

  l = nm = 0;

	if (m < n)
		nrerror("SVDCMP: You must augment A with extra zero rows");
	init_svdcmp(n, 1, rbuf);
	/* Householder reduction to bidiagonal form */
	for (i = 1; i <= n; i++) {
		l = i + 1;
		(glob.rv1)[i] = scale * g;
		g = s = scale = 0.0;
		if (i <= m) {
			for (k = i; k <= m; k++) 
				scale += fabs(a[k][i]);
			if (scale) {
				for (k = i; k <= m; k++) {
					a[k][i] /= scale;
					s += a[k][i] * a[k][i];
				}
				f = a[i][i];
				g = -SIGN(sqrt(s), f);
				h = f * g - s;
				a[i][i] = f - g;
				if (i != n) {
					for (j = l; j <= n; j++) {
						for (s = 0.0, k = i; k
						    <= m; k++) 
							s += a[k][i] *
							    a[k][j];
						f = s / h;
						for (k = i; k <= m; k++) 
							a[k][j] += f *
							    a[k][i];
					}
				}
				for (k = i; k <= m; k++) 
					a[k][i] *= scale;
			}
		}
		w[i] = scale * g;
		g = s = scale = 0.0;
		if (i <= m && i != n ) {
			for (k = l; k <= n; k++) 
				scale += fabs(a[i][k]);
			if (scale) {
				for (k = l; k <= n; k++) {
					a[i][k] /= scale;
					s += a[i][k] * a[i][k];
				}
				f = a[i][l];
				g =  -SIGN(sqrt(s), f);
				h = f * g - s;
				a[i][l] = f - g;
				for (k = l; k <= n; k++) 
					(glob.rv1)[k] = a[i][k] / h;
				if (i != m) {
					for (j = l; j <= m; j++) {
						for (s = 0.0, k = l; k
						    <= n; k++)
							s += a[j][k] *
							    a[i][k];
						for (k = l; k <= n; k++)
							a[j][k] += s *
							    (glob.rv1)[k];
					}
				}
				for (k = l; k <= n; k++) 
					a[i][k] *= scale;
			}
		}
		anorm = MAX(anorm, (fabs(w[i]) + fabs((glob.rv1)[i])));
	}

	if (0.0 + anorm != anorm) /* Test if NaN */
	  nrerror("Nonsense matrix elements in svdcmp!");

	/* Accumulation of right-hand transformations */
	for (i = n; i >= 1; i--) {
		if (i < n) {
			if (g) {
				for (j = l; j <= n; j++) 
                              /* Double division to avoid possible underflow */
					v[j][i] = (a[i][j] / a[i][l]) /
					    g;
				for (j = l; j <= n; j++) {
					for (s = 0.0, k = l; k <= n; k++) 
						s += a[i][k] * v[k][j];
					for (k = l; k <= n; k++) 
						v[k][j] += s * v[k][i];
				}
			}
			for (j = l; j <= n; j++) 
				v[i][j] = v[j][i] = 0.0;
		}
		v[i][i] = 1.0;
		g = (glob.rv1)[i];
		l = i;
	}
	/* Accumulation of left-hand transformations */
	for (i = n; i >= 1; i--) {
		l = i + 1;
		g = w[i];
		if (i < n)
			for (j = l; j <= n; j++) 
				a[i][j] = 0.0;
		if (g) {
			g = 1.0 / g;
			if (i != n) {
				for (j = l; j <= n; j++) {
					for (s = 0.0, k = l; k <= m; k++) 
						s += a[k][i] * a[k][j];
					f = (s / a[i][i]) * g;
					for (k = i; k <= m; k++) 
						a[k][j] += f * a[k][i];
				}
			}
			for (j = i; j <= m; j++) 
				a[j][i] *= g;
		} else {
			for (j = i; j <= m; j++) 
				a[j][i] = 0.0;
		}
		++a[i][i];
	}
	/* Diagonalization of the bidiagonal form */
	for ( k = n; k >= 1; k--) {    /* Loop over singular values */
		for (its = 1; its <= ITERNUM; its++) {
			/* Loop over allowed iterations */
			flag = 1;
			for (l = k; l >= 1; l--) {     /* Test for splitting */
			        nm = l - 1;            /* glob.rv1[1] = 0! */
				if (fabs((glob.rv1)[l]) + anorm == anorm) {
					flag = 0;
					break;
				}
				if (fabs(w[nm]) + anorm == anorm) 
					break;
			}
			if (flag) {
				c = 0.0;
				s = 1.0;
				for (i = l; i <= k; i++) {
					f = s * (glob.rv1)[i];
					if (fabs(f) + anorm != anorm) {
						g = w[i];
						h = PYTHAG(f, g);
						w[i] = h;
						h = 1.0 / h;
						c = g * h;
						s = (-f * h);
						for (j = 1; j <= m; j++) {
							y = a[j][nm];
							z = a[j][i];
							a[j][nm] = y *
							    c + z * s;
							a[j][i] =  z *
							    c - y * s;
						}
					}
				}
			}
			z = w[k];
			if (l == k) {			   /* Convergence */
				if (z < 0.0) {   /* Sing. value is made nonnegtive */
					w[k] = -z;
					for (j = 1; j <= n; j++) 
						v[j][k] = (-v[j][k]);
				}
				break;
			}
			if (its == ITERNUM)
			{ nrerror("No convergence in ITERNUM SVDCMP iterations");
                          return(0);
			  }
			x = w[l];
			nm = k - 1;
			y = w[nm];
			g = (glob.rv1)[nm];
			h = (glob.rv1)[k];
			/*	      f= ((y-z)*(y+z) + (g-h)*(g+h))/(2.0*h*y);
*/
			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0
			    *h);
			f /= y;
			g = PYTHAG(f, 1.0);
			f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g,
			     f))) - h)) / x;
			/* Next QR transformation */
			c = s = 1.0;
			for (j = l; j <= nm; j++) {
				i = j + 1;
				g = (glob.rv1)[i];
				y = w[i];
				h = s * g;
				g = c * g;
				z = PYTHAG(f, h);
				(glob.rv1)[j] = z;
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = g * c - x * s;
				h = y * s;
				y = y * c;
				for (jj = 1; jj <= n; jj++) {
					x = v[jj][j];
					z = v[jj][i];
					v[jj][j] = x * c + z * s;
					v[jj][i] = z * c - x * s;
				}
				z = PYTHAG(f, h);
				w[j] = z; /* Rotation can be arbitrary if z=0 */
				if (z>MINDOUBLE) {
					z = 1.0 / z;
					c = f * z;
					s = h * z;
				}
				f = (c * g) + (s * y);
				x = (c * y) - (s * g);
				for (jj = 1; jj <= m; jj++) {
					y = a[jj][j];
					z = a[jj][i];
					a[jj][j] = y * c + z * s;
					a[jj][i] = z * c - y * s;
				}
			}
			(glob.rv1)[l] = 0.0;
			(glob.rv1)[k] = f;
			w[k] = x;
		}
	}
  return(1);
}




