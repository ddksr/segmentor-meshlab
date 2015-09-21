// MATRIX.C

#include "matrix.h"

#include <iostream>
#include <math.h>
#include <stdio.h>

// Eigenvalue Decomposition routine

/* Function for computing eigenvalues and eigenvectors, no sorting

            call: eigendcmp (a, n, d, v, &nr);

            Input: a : pointer to the input array of type double
                   n : input (output) array size: n x n (int)

            Output:
                   d : pointer to the output array holding eigenvalues (double)
                   v : pointer to the output array holding eigenvectors (double)
                       v[i][j] column eigenvectors j, j = 0, 1, ..., n-1    
                  nr : number of iterations performed (int).

*/

int eigendcmp( double *a, int n, double *d, double *v, int *nr)
{

int   nrot;

int   j, iq, ip, i;

double tresh,  theta,  tau,  t,  sm,  s,  h,  g,  c;

double b[256], z[256];


  for (ip = 0; ip < n; ip++){

    for (iq = 0; iq < n; iq++){

      v[ip * n + iq] = 0.0;
    }
    v[ip * n  + ip] = 1.0;
  }

  for (ip = 0; ip < n; ip++){

    b[ip] = a[ip * n + ip];
    d[ip] = b[ip];
    z[ip] = 0.0;
  }

  nrot = 0;

  for (i = 1; i < 50; i++){

    sm = 0.0;

    for ( ip = 0; ip < n-1; ip++){

      for ( iq = ip+1; iq <  n; iq++){

        sm += fabs( a[ip * n + iq]);

      }

    }

    if (sm == 0.0){
       *nr = nrot;
       return 0;
    }

    if (i < 6){

      tresh = (double)0.2 * sm / (double)(n * n);
    }

    else{

      tresh = 0.0;

    }

    for (ip = 0; ip < n-1; ip++){

      for (iq = ip+1; iq < n; iq++){

        g = 100.0 * fabs(a[ip * n + iq]);

        if ((i > 6) && ((fabs(d[ip]) + g) == fabs(d[ip]))
                    && ((fabs(d[iq]) + g) == fabs(d[iq]))){
           a[ip * n + iq] = 0.0;
        }

        else if (fabs(a[ip * n + iq]) > tresh){

          h = d[iq] - d[ip];

          if ((fabs(h) + g) == fabs(h)){

            t = a[ip * n + iq] / h;
	  }

          else{
            theta = 0.5 * h / a[ip * n + iq];
            t     = 1.0 /(fabs(theta) + sqrt(1.0 + theta * theta));
            if (theta < 0.0)  t = -t;
	  }

          c   = 1.0 / sqrt(1 + t * t);
          s   = t * c;
          tau = s / (1.0 + c);
          h   = t * a[ip * n + iq];

          z[ip] = z[ip] - h;
          z[iq] = z[iq] + h;
          d[ip] = d[ip] - h;
          d[iq] = d[iq] + h;

          a[ip * n + iq] = 0.0;

          for (j = 0; j <= ip-1; j++){

            g = a[j * n + ip];
            h = a[j * n + iq];
            a[j * n + ip] = g - s * (h + g * tau);
            a[j * n + iq] = h + s * (g - h * tau);
	  }

          for (j = ip+1; j <= iq-1; j++){

            g = a[ip * n + j];
            h = a[j  * n + iq];
            a[ip * n + j ] = g - s * (h + g * tau);
            a[j  * n + iq] = h + s * (g - h * tau);
	  }

          for (j = iq+1; j < n; j++){
            g = a[ip * n + j];
            h = a[iq * n + j];
            a[ip * n + j] = g - s * (h + g * tau);
            a[iq * n + j] = h + s * (g - h * tau);
	  }

          for (j = 0; j < n; j++){
            g = v[j * n + ip];
            h = v[j * n + iq];
            v[j * n + ip] = g - s * (h + g * tau);
            v[j * n + iq] = h + s * (g - h * tau);
	  }
          nrot++;

	}
      }
    }

    for (ip = 0; ip < n; ip++){

      b[ip] = b[ip] + z[ip];
      d[ip] = b[ip];
      z[ip] = 0.0;
    }
  }
  printf ("routine JACOBI: 50 iterations should not happen\n");
  *nr = nrot;
  return -1;

} /* end eigen */


matrix::matrix(int h1, int w1)
{ w = w1;
  h = h1;
  elem = new double[w*h];
  //  std::cout << "Cons ";
  }

matrix::matrix(const matrix &a)
{ int i,n;

  w = a.w; h = a.h;
  elem = new double[n = w * h];
  for (i = 0; i < n; i++) elem[i] = a.elem[i];
  // std::cout << "Cons ";
  }

matrix& matrix::operator=(const matrix& a)
{ int i,n;

  if (this != &a)
  { 
    delete [] elem;
    w = a.w; h = a.h;
    elem = new double[n = w * h];
    for (i = 0; i < n; i++) elem[i] = a.elem[i];
    }
    return(*this);
  }

matrix& matrix::operator+=(matrix a)
{ int i,n;

  if (w == a.w && h == a.h)
  { n = w * h;
    for (i = 0; i < n; i++) elem[i]+=a.elem[i];
    } else
  { std::cout << "Not compatible types of matrix in operator +";
    }
  return (*this);
  }

matrix& matrix::operator-=(matrix a)
{ int i,n;

  if (w == a.w && h == a.h)
  { n = w * h;
    for (i = 0; i < n; i++) elem[i]-=a.elem[i];
    } else
  { std::cout << "Not compatible types of matrix in operator -";
    }
  return (*this);
  }

matrix& matrix::operator/=(double x)
{ int i,n;
  
  n = w * h;
  if (x != 0.0)
    for (i = 0; i < n; i++) elem[i] /= x;
  else
    std::cout << "Division by zero in matrix operator /=";
  return(*this);
  }

matrix& matrix::operator*=(matrix a)
{ double *d;
  int i,j,k,n;
  if (w == a.h)
  { d = new double[h * a.w];
    for (i = 0; i < a.w; i++)
      for (j = 0; j < h; j++)
      { n = i + j * a.w;
	d[n] = 0;
	for (k = 0; k < w; k++)
	  d[n] = d[n] + el(j,k) * a.el(k,i);
	}
    w = a.w;
    delete [] elem;
    elem = d;
    } else
  { std::cout << "Not compatible types of matrix in operator *";
    }
  return(*this);
  }

matrix matrix::operator-()
{ matrix r = *this;
  int i,n = w * h;
  for (i = 0; i < n; i++) r.elem[i] = -r.elem[i];
  return(r);
  }

matrix matrix::transpose()
{ matrix r(w,h);
  int i,j;

  for (i = 0; i < w; i++)
    for (j = 0; j < h; j++)
      r.el(i,j) = el(j,i);
  return(r);
  }


void matrix::exch_rows(int n1, int n2)
{ int i;
  double t;

  for (i = 0; i < w; i++)
  { t = el(n1,i);
    el(n1,i) = el(n2,i);
    el(n2,i) = t;
    }
  }

void matrix::exch_cols(int n1, int n2)
{ int j;
  double t;

  for (j = 0; j < h; j++)
  { t = el(j,n1);
    el(j,n1) = el(j,n2);
    el(j,n2) = t;
    }
  }

void matrix::multiply_row(int n, double d)
{ int i;
  for (i = 0; i < w; i++) el(n,i)*=d;
  }

void matrix::multiply_col(int n, double d)
{ int j;
  for (j = 0; j < h; j++) el(j,n)*=d;
  }

void matrix::add_mul_row(int n1, int n2, double d)
{ int i;
  for (i = 0; i < w; i++) el(n2,i)+=d*el(n1,i);
  }

void matrix::add_mul_col(int n1, int n2, double d)
{ int j;
  for (j = 0; j < w; j++) el(j,n2)+=d*el(j,n1);
  }

void matrix::set(int y, int x, double d)
{ el(y,x) = d;
  }

matrix& matrix::numeric_round()
{ int i,n = w*h;

  for (i = 0; i < n; i++)
  { if (fabs(elem[i]-floor(elem[i])) < 1e-6) elem[i] = floor(elem[i]);
      else if (fabs(elem[i]-ceil(elem[i])) < 1e-6) elem[i] = ceil(elem[i]);
    }
  return(*this);
  }

void matrix::set_submatrix(matrix& a, int r, int c)
{ int i,j,dw,dh;

  dw = r + a.nrows(); 
  dh = c + a.ncols();

  for (i = r; i < dw; i++)
    for (j = c; j < dh; j++)
      el(i,j) = a.el(i-r,j-c);
  } 

void matrix::load(FILE *f)
{ fread(&w,sizeof(int),1,f);
  fread(&h,sizeof(int),1,f);
  delete [] elem;
  elem = new double[w*h];
  fread(elem,sizeof(double),w*h,f);
  }  
  
void matrix::save(FILE *f)
{ fwrite(&w,sizeof(int),1,f);
  fwrite(&h,sizeof(int),1,f);
  fwrite(elem,sizeof(double),w*h,f);
  }  
  
matrix::~matrix()
{ delete [] elem;
  // std::cout << "Dest ";
  }

matrix operator+(matrix a, matrix b)
{ matrix r = a;
  r += b;
  return(r);
  }

matrix operator-(matrix a, matrix b)
{ matrix r = a;
  r -= b;
  return(r);
  }

matrix operator*(matrix a, matrix b)
{ matrix r = a;
  r *= b;
  return(r);
  }

matrix operator*(matrix a, double d)
{ int i,n = a.w * a.h;
  matrix r = a;
  for (i = 0; i < n; i++) r.elem[i]*=d;
  return(r);
  }

matrix operator*(double d,matrix a)
{ int i,n = a.w * a.h;
  matrix r = a;
  for (i = 0; i < n; i++) r.elem[i]*=d;
  return(r);
  }

std::ostream& operator<<(std::ostream& s, const matrix a)
{ int i,j;

  for (j = 0; j < a.nrows(); j++)
  { for (i = 0; i < a.ncols(); i++)
      s << a.el(j,i) << ' ';
    s << '\n';
    }
  return(s);
  }

symatrix::symatrix(int n):matrix(n,n)
{
  }

void symatrix::identity()
{ int i,j;
  for (j = 0; j < h; j++)
    for (i = 0; i < w; i++)
      if (i == j) el(i,j) = 1;
	else el(i,j) = 0;
  }


void symatrix::solve(vector &x, vector &b)
{ double *aa, *vv, *ww;
   symatrix V(w), W(w),F(w);
  int i,j,k;

  aa = new double[w*h];
  vv = new double[w*h];
  ww = new double[w];
  
  for (i = 0; i < w; i++)
     for (j = 0; j < h; j++)
        aa[i*w+j] = el(i,j);
        
  eigendcmp(aa,w,ww,vv,&k);
  
  for (i = 0; i < w; i++)
  { for (j = 0; j < h; j++)
     { V.el(i,j) = vv[i*w+j];
        W.el(i,j) = 0.0;
        }
     W.el(i,i) = ww[i];
     }

//  std::cout << "\n\nV = " << V << "\n\nW = " << W << "\n\n";
//  std::cout << "V*V' = " << V * V.transpose() << "\n\n";

   for (i = 0; i < w; i++) 
      if (fabs(W.el(i,i)) > 0.0001) W.el(i,i) = 1.0 / W.el(i,i);

  F = V * W * V.transpose();
  x = F * b;


//  std::cout << "F = " << F << "\n\n";
//  std::cout << "A = " << (*this) << "\n\ninv(A) = " << F << "\n\n";
//  std::cout << "A*inv(A) = " << (*this) * F << "\n\n";

  delete [] aa;
  delete [] vv;
  delete [] ww;
  }


symatrix symatrix::inverse(int& sing)
{ symatrix inv = *this, tm = *this, ch = *this;
  int i,j;
  double d;
  inv.identity();

  sing = 0;
  for (i = 0; i < tm.w && !sing; i++)
  { j = i;
    while (j < tm.h && tm.el(j,i) == 0) j++;
    if (j >= tm.h)
      { std::cout << "\nWarning:Singular matrix" << '\n' << *this;
        sing = 1;
        }
      else
      { tm.exch_rows(i,j);
        inv.exch_rows(i,j);           // rows or cols this is the question?
	for (j = i + 1; j < tm.h; j++)
	{ d = -tm.el(j,i)/tm.el(i,i);
	  tm.add_mul_row(i,j,d);
	  inv.add_mul_row(i,j,d);
	  }
	d = 1 / tm.el(i,i);
	tm.multiply_row(i,d);
	inv.multiply_row(i,d);
	}
    }
  if (!sing)
    for (i = tm.w - 1; i > 0; i--)
      for (j = 0; j < i; j++)
      { d = -tm.el(j,i)/tm.el(i,i);
	tm.add_mul_row(i,j,d);
	inv.add_mul_row(i,j,d);
	}

  ch = *this * inv;
  ch.numeric_round();

  for (i = 0; i < tm.w; i++)
    for (j = 0; j < tm.h; j++)
    { if (i == j) d = 1.0;
        else d = 0.0;
      if (ch.el(i,j) != d) sing = 2;
      }
  
//  if (sing == 2) std::cout << "Inverse not computed correctly: \n" << ch << "\n" << *this << "\n" << inv;

  return(inv);
  }

double symatrix::determinant()
{ symatrix tm = *this;
  double det = 1.0,d;
  int i,j,sing = 0;

  for (i = 0; i < tm.w && !sing; i++)
  { j = i;
    while (j < tm.h && tm.el(j,i) == 0) j++;
    if (j >= tm.h)
    { sing = 1;
      }
      else
      { tm.exch_rows(i,j);
        if (i != j) det = -det;
	for (j = i + 1; j < tm.h; j++)
	{ d = -tm.el(j,i)/tm.el(i,i);
	  tm.add_mul_row(i,j,d);
	  }
	d = tm.el(i,i);
        det *= d;
        }
    }
  if (sing) det = 0.0;
  return(det);
  }

symatrix::symatrix(const symatrix& a):matrix(a.nrows(),a.nrows())
{ int i,j,x = a.ncols(),y = a.nrows();

  for (j = 0; j < y; j++)
    for (i = 0; i < x; i++) el(j,i) = a.el(j,i);
 
  }

symatrix symatrix::operator=(matrix a)
{ int i,j,x = a.ncols(),y = a.nrows();

  if (x != y)
  { std::cout << "\nMatrix attempted to assign to symatrix!";
    }
  else
  if (this != &a)
    { delete [] elem;
      h = y; w = x;
      elem = new double[h*w];
      for (j = 0; j < h; j++)
	for (i = 0; i < w; i++) el(j,i) = a.el(j,i);
      }

  return(*this);
  }

vector::vector(int n):matrix(n,1)
{
  }

vector::vector(const vector& a):matrix(a.nrows(),1)
{ int i,n = a.nrows();

  for (i = 0; i < n; i++) el(i) = a.el(i);
  }

vector vector::operator=(matrix a)
{ int i,n = a.nrows();

  if (a.ncols() != 1)
  { std::cout << "\nMatrix attempted to assign to vector";
    } else
  if (this != &a)
    { delete [] elem;
      h = n; w= 1;
      elem = new double[h];
      for (i = 0; i < h; i++) el(i) = a.el(i,0);
      }
  return(*this);
  }

double vector::norm()
{ return(sqrt(this->inner(*this)));
  }

vector vector::normalize()
{ double normal = norm();
  vector v = *this;

  if (normal > 0) v.multiply_col(0,1/normal);

  return(v);
  }

double& vector::el(int n) const
{ return(matrix::el(n,0));
  }  
  
double operator*(vector a, vector b)
{ double sum = 0;
  int i,n = a.nrows();
  if (n == b.nrows())
  { for (i = 0; i < n; i++) sum += a.el(i) * b.el(i);
    } else
  { std::cout << "\nIncompatible types of vectors in operator *";
    }
  return(sum);
  }

double vector::inner(vector a)
{ double sum = 0.0;
  int i,n = nrows();
  
  if (n == a.nrows())
  { for (i = 0; i < n; i++) sum += el(i) * a.el(i); 
    } else
  { std::cout << "\nIncompatible types of vectors in inner";
    }
  return(sum);
  }

double vector::distance(vector v)
{ int i,n = nrows();
  double sum = 0.0;

  if (n == v.nrows())
    for (i = 0; i < n; i++) sum += (el(i) - v.el(i)) * (el(i) - v.el(i));
      else 
	std::cout << "\nIncompatible types of vectors in distance";

  return(sqrt(sum));
  }

vector3D::vector3D():vector(3)
{
  }

vector3D::vector3D(const vector3D& a):vector(3)
{ int i;

  for (i = 0; i < 3; i++) el(i) = a.el(i);
  }

vector3D& vector3D::operator=(matrix a)
{ int i,x = a.ncols(),y = a.nrows();
  if (x != 1 || y != 3)
  { std::cout << "\nMatrix attempted to assign to vector3D";
    } else
    if (this != &a)
    { 
      delete [] elem;
      h = 3; w = 1;
      elem = new double[3];
      for (i = 0; i < 3; i++) elem[i] = a.el(i,0);
      }
  return(*this);
  }

vector3D vector3D::outer(vector3D& a)
{ vector3D v;

  v.el(0) = el(1)*a.el(2) - el(2)*a.el(1);
  v.el(1) = el(2)*a.el(0) - el(0)*a.el(2);
  v.el(2) = el(0)*a.el(1) - el(1)*a.el(0);

  return(v);
  }

hmatrix::hmatrix() : symatrix(4)
{ identity();
  }
  
hmatrix::hmatrix(const hmatrix& a) : symatrix(a)
{
  }
  
hmatrix hmatrix::operator=(matrix a)
{ int i,j,x = a.ncols(), y = a.nrows();
  if (x != 4 || y != 4)
  { std::cout << "\nMatrix attempted to assign to hmatrix";
    } else
  if (this!=&a)
  { for (i = 0; i < 4; i++)
      for (j = 0; j < 4; j++) el(i,j) = a.el(i,j);
    }
  return(*this);
  }  
  
  
// written by Ales Jaklic
// translation relative to the relative frame

void hmatrix::translate_r(double x, double y, double z) {
  
  hmatrix tmp;
  
  tmp.el(0,3) = x;
  tmp.el(1,3) = y;
  tmp.el(2,3) = z;
  
  *this = *this * tmp;
}

// translation relative to the global frame

void hmatrix::translate_g(double x, double y, double z) {
  
  hmatrix tmp;
  
  tmp.el(0,3) = x;
  tmp.el(1,3) = y;
  tmp.el(2,3) = z;
  
  *this = tmp * (*this);
}

// rotation about the x axis of the relative frame

void hmatrix::rotate_xr(double alpha) {
  
  hmatrix tmp;
  
  tmp.el(1,1) = cos(alpha);
  tmp.el(1,2) = -sin(alpha);
  tmp.el(2,1) = sin(alpha);
  tmp.el(2,2) = cos(alpha);
  
  
  *this = (*this) * tmp;
} 

// rotation about the x axis of the global frame

void hmatrix::rotate_xg(double alpha) {
  
  hmatrix tmp;
  
  tmp.el(1,1) = cos(alpha);
  tmp.el(1,2) = -sin(alpha);
  tmp.el(2,1) = sin(alpha);
  tmp.el(2,2) = cos(alpha);
  
  
  *this = tmp * (*this);
} 	                     

// rotation about the y axis of the relative frame

void hmatrix::rotate_yr(double alpha) {
  
  hmatrix tmp;
  
  tmp.el(0,0) = cos(alpha);
  tmp.el(0,2) = sin(alpha);
  tmp.el(2,0) = -sin(alpha);
  tmp.el(2,2) = cos(alpha);
  
  
  *this = (*this) * tmp;
}

// rotation about the y axis of the general frame

void hmatrix::rotate_yg(double alpha) {
  
  hmatrix tmp;
  
  tmp.el(0,0) = cos(alpha);
  tmp.el(0,2) = sin(alpha);
  tmp.el(2,0) = -sin(alpha);
  tmp.el(2,2) = cos(alpha);
  
  
  *this = tmp * (*this);
}

// rotation about the z axis of the relative frame

void hmatrix::rotate_zr(double alpha) {
  
  hmatrix tmp;
  
  tmp.el(0,0) = cos(alpha);
  tmp.el(0,1) = -sin(alpha);
  tmp.el(1,0) = sin(alpha);
  tmp.el(1,1) = cos(alpha);
  
  
  *this = (*this) * tmp;
}

// rotation about the z axis of the general frame

void hmatrix::rotate_zg(double alpha) {
  
  hmatrix tmp;
  
  tmp.el(0,0) = cos(alpha);
  tmp.el(0,1) = -sin(alpha);
  tmp.el(1,0) = sin(alpha);
  tmp.el(1,1) = cos(alpha);
  
  
  *this = tmp * (*this);
}
  
  
/* How to use classes matrix, symatrix, vector and vector3D

void main()
{
  matrix A(2,3);
  matrix J(3,1);
  symatrix D(4);
  vector C(2);
  symatrix S1(4);
  A.el(0,0) = A.el(1,1) = 1;
  A.el(1,0) = 2;
  A.el(0,1) = -2; A.el(0,2) = -2; A.el(1,2) = 3;
  J.el(0,0) = 1; J.el(1,0) = 0;
  J.el(2,0) = 5;
  matrix B = A.transpose();
  std::cout << '+';
  C.el(0) = 1.5;
  C.el(1) = 4;
  S1.identity();
  S1.multiply_row(3,3);
  S1.exch_rows(1,3);
  S1.add_mul_row(3,1,2);
  std::cout << "\n\n\n\n\n\n" << C;
  C = A * J;
  std::cout << '\n' << C;
//  std::cout << -A << -B << A*C << '\n' << C*C << D << C.normalize() << '\n';
//  std::cout << C.norm() << '\n' << S1 << S1.inverse() << S1*S1.inverse();

  D.el(0,0) = 2; D.el(0,1) = 1; D.el(0,2) = -3; D.el(0,3) = 0;
  D.el(1,0) = 1; D.el(1,1) = 1; D.el(1,2) = 3; D.el(1,3) = -2;
  D.el(2,0) = -1; D.el(2,1) = 0; D.el(2,2) = 4; D.el(2,3) = 1;
  D.el(3,0) = 1; D.el(3,1) = 1; D.el(3,2) = 7; D.el(3,3) = -1;

//  std::cout << D << '\n' << D.inverse() << '\n' << (D*D.inverse()).numeric_round() <<
//  '\n' << (D.inverse()*D).numeric_round();
  } */
