// matrix.H

#ifndef LIBSEGMENTOR_MATRIX
#define LIBSEGMENTOR_MATRIX 1

#include <iostream>
#include <stdio.h>

class vector;

class matrix
{
protected:
  int w,h;
  double *elem;
public:
  matrix(int h1, int w1);      // constructor;
  matrix(const matrix &a);     // copy constructor;
  matrix& operator=(const matrix& a);
  int is_in(int x, int y) const { return(x >= 0 && x < w && y >= 0 && y < h); }
  double& el(int y, int x) const 
  { if (is_in(x,y)) return(elem[x+w*y]); 
      else
	{ std::cout << "\nOut of range in matrix::el, element " << x << ',' << y << " Matrix(" << w << "," << h << ")";
      return(elem[0]);
      }
    }
  matrix& operator+=(matrix a);
  matrix& operator-=(matrix a);
  matrix& operator*=(matrix a);
  matrix& operator/=(double x);
  matrix operator-();
  matrix transpose();
  inline int nrows() const { return(h); }
  inline int ncols() const { return(w); }
  void exch_rows(int n1, int n2);
  void exch_cols(int n1, int n2);
  void multiply_row(int n, double d);
  void multiply_col(int n, double d);
  void add_mul_row(int n1, int n2, double d);
    //add row n1 multiplied by d to row n2
  void add_mul_col(int n1, int n2, double d);
    //add col n1 multiplied by d to col n2
  void set(int y, int x, double d);
  matrix& numeric_round();
  friend matrix operator*(matrix a, double d);
  friend matrix operator*(double d, matrix a);
  matrix& makenull();
  void load(FILE *f);
  void save(FILE *f);
  void set_submatrix(matrix& a, int r, int c); 
  virtual ~matrix();
  };

matrix operator+(matrix a, matrix b);
matrix operator-(matrix a, matrix b);
matrix operator*(matrix a, matrix b);
std::ostream& operator<<(std::ostream& s, const matrix a);

class symatrix : public matrix
{
public:
  symatrix(int n);
  void identity();
  symatrix inverse(int& sing);  
  symatrix(const symatrix& a);
  symatrix operator=(matrix a);
  double determinant(); 
  void solve(vector &x, vector &b);         // solves Ax=b, A == this 
  };

class vector : public matrix
{
public:
  vector(int n);
  vector(const vector& a);
  vector operator=(matrix a);
  double norm();
  vector normalize();
  double& el(int n) const;
  double inner(vector v);
  double distance(vector v);
  };

double operator*(vector a, vector b);

class vector3D : public vector
{
public:
  vector3D();
  vector3D(const vector3D& a);
  vector3D& operator=(matrix a);
  vector3D outer(vector3D& a);
  };
  
class hmatrix : public symatrix
{ 
public:
  hmatrix();
  hmatrix(const hmatrix& a);
  hmatrix operator=(matrix a);
  void translate_g(double, double, double);
  void translate_r(double, double, double);
  void rotate_xr(double);
  void rotate_xg(double);
  void rotate_yr(double);
  void rotate_yg(double);
  void rotate_zr(double);
  void rotate_zg(double);
  void unit() { identity(); }
  }; 
  
#endif








