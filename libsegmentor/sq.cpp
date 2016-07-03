// sq.C

#include "sq.h"
#include "state.h"

#include <math.h>
#include <iostream>

extern "C" {
#include "sq_recover.h"
 double sq_mindist(double a1, double a2, double a3, double e1, double e2, double x0, double y0, double z0);
}

double spow(double x, double y) 
{
  if (x < 0.0)
    return(-pow(-x, y));
  else if (x == 0.0)
    return(0);
  else
    return(pow(x, y));
 
}

int recover_params(sq* model, vect* list, int no) {
  if (State::lookup->on && State::lookup->isSet) {
	// Reset SQ params to those from lookup
	model->a1 = State::lookup->a;
	model->a2 = State::lookup->b;
	model->a3 = State::lookup->c;
	model->e1 = State::lookup->e1;
	model->e2 = State::lookup->e2;

	double _kx = 0, _ky = 0;
	return recover_search(list, no,
						  &model->a1, &model->a2, &model->a3, &model->e1, &model->e2,
						  &model->px, &model->py, &model->pz, &model->phi, &model->theta, &model->psi,
						  &_kx, &_ky, RECOVER_SQ);
  }
  return recover(list, no,
				 &model->a1, &model->a2, &model->a3, &model->e1, &model->e2,
				 &model->px, &model->py, &model->pz, &model->phi, &model->theta, &model->psi);
}

sq::sq(region &r) {

  #define MAX_LEN 25000
  
  struct vect list[MAX_LEN];
  int x, y, min_x, min_y, max_x, max_y, no;
  struct point p;  
  int sing;

  no = 0;

  min_x = r.minx;
  min_y = r.miny;
  max_x = r.maxx;
  max_y = r.maxy;

  for(x = min_x; x <= max_x; x++)
    for(y = min_y; y <= max_y; y++)
      if (r.included(x, y))
	if (no < MAX_LEN) 
	{ p = r.get_point(x,y);
          list[no].x = p.x;
          list[no].y = p.y;
          list[no++].z = p.z;
          }
	else
	  std::cout << "Over the limit of SQ points";

  double t[4][4];
 
  estimate(list, no, t);
  convert(list, no, t, &a1, &a2, &a3, &e1, &e2, &px, &py, &pz, &phi, &theta, &psi);
  //std::cout << "Estimate: \n";
  //print();
  //std::cout << "------------\n";
  
  printf("\nSQ recovering from scratch    %d points", no); fflush(stdout);
  if (!recover_params(this, list, no))
  { printf("\nSQ recover procedure returned error!"); fflush(stdout);
    }

  printf("\nSQ recover ended"); fflush(stdout);
  g_from_l.translate_g(px, py, pz);
  g_from_l.rotate_zr(phi);
  g_from_l.rotate_yr(theta);
  g_from_l.rotate_zr(psi);
  
  l_from_g = g_from_l.inverse(sing);
}

// creates an SQ model based on region and given initial estimate

sq::sq(sq *m, region& r) {


  #define MAX_LEN 25000

  struct vect list[MAX_LEN];
  int x, y, min_x, min_y, max_x, max_y, no;
  struct point p;
  int sing;
  
  a1 = m->a1; a2 = m->a2; a3 = m->a3; e1 = m->e1; e2 = m->e2;
  phi = m->phi; theta = m->theta; psi = m->psi; px = m->px; py = m->py; 
  pz = m->pz;
  
  no = 0;

  min_x = r.minx;
  min_y = r.miny;
  max_x = r.maxx;
  max_y = r.maxy;

  for(x = min_x; x <= max_x; x++)
    for(y = min_y; y <= max_y; y++)
      if (r.included(x, y))
	if (no < MAX_LEN) 
        { p = r.get_point(x,y);
          list[no].x = p.x;
          list[no].y = p.y;
          list[no++].z = p.z;
          }
	else
	  std::cout << "Over the limit of SQ points";

   
 
  printf("\nSQ recovering from old parameters"); fflush(stdout);
  if (!recover_params(this, list, no))
  { printf("\nSQ recover procedure returned error!"); fflush(stdout);
    }

  printf("\nSQ recover ended"); fflush(stdout);
  g_from_l.translate_g(px, py, pz);
  g_from_l.rotate_zr(phi);
  g_from_l.rotate_yr(theta);
  g_from_l.rotate_zr(psi);
  
  l_from_g = g_from_l.inverse(sing);
}

sq::sq(FILE *f)
{ fread(&a1, 1, sizeof(double), f);
  fread(&a2, 1, sizeof(double), f);
  fread(&a3, 1, sizeof(double), f);
  fread(&e1, 1, sizeof(double), f);
  fread(&e2, 1, sizeof(double), f);
  fread(&phi, 1, sizeof(double), f);
  fread(&theta, 1, sizeof(double), f);
  fread(&psi, 1, sizeof(double), f);
  fread(&px, 1, sizeof(double), f);
  fread(&py, 1, sizeof(double), f);
  fread(&pz, 1, sizeof(double), f);
  l_from_g.load(f);
  g_from_l.load(f);
  }

double sq::map_eta(double eta) {

  double x = a2 * spow(cos(eta), 1/e1);
  double y = a3 * spow(sin(eta), 1/e1);
  
  return (atan2(y, x));
}

double sq::map_omega(double omega) {
  double x = a2 * spow(cos(omega), 1/e2);
  double y = a1 * spow(sin(omega), 1/e2);

  return (atan2(y, x));
}

vector sq::normal(double eta, double omega) const {
 
  vector result(4);

  result.el(0) = spow(cos(eta), 2.0 - e1) * spow(cos(omega), 2.0 - e2) / a1;
  result.el(1) = spow(cos(eta), 2.0 - e1) * spow(sin(omega), 2.0 - e2) / a2;
  result.el(2) = spow(sin(eta), 2.0 - e1) / a3;
  result.el(3) = 1.0;
  
  result = (1.0 / result.norm()) * (matrix)result;


  return result;
}

vector sq::r(double eta, double omega) const {
 
  vector result(4);

  result.el(0) = a1 * spow(cos(eta), e1) * spow(cos(omega), e2);
  result.el(1) = a2 * spow(cos(eta), e1) * spow(sin(omega), e2);
  result.el(2) = a3 * spow(sin(eta), e1);
  result.el(3) = 1.0;
  
  return result;
}

double sq::abs_signed_distance(struct point& p)
{ vector tmp(4),r(4);
  double delta;
  tmp.el(0) = p.x;
  tmp.el(1) = p.y;
  tmp.el(2) = p.z;
  tmp.el(3) = 1.0;
  
  r = l_from_g * tmp;
  
  delta = sqrt(r * r) * ((1.0 - pow(f(r.el(0),r.el(1),r.el(2)),-0.5)));
  
  return(delta);
  }
  
double sq::f(double x, double y, double z) const {
  

  // thid double pow is to prevent negative arguments to pow when the exponent is negative
  // do not try to optimize

  double a = pow(pow((x / a1), 2), 1/e2);
  double b = pow(pow((y / a2), 2), 1/e2);   
  double c = pow(pow((z / a3), 2), 1/e1);
  
  return pow(pow(a + b, e2/e1) + c, e1);
}

void sq::sq_draw(int no_hidding) {

  #define PI 3.14159265359
  #define N

  int i, j;  
  double omega, eta;
  vector new_r(4),old_r(4),new_r1(4),old_r1(4);

  new_r.el(0) = new_r.el(1) = new_r.el(2) = old_r.el(0) = old_r.el(1) = old_r.el(2) = 0;
  new_r.el(3) = old_r.el(3) = 1;
  
  vector view_g(4);
  view_g.el(0) = 0;
  view_g.el(1) = 0;
  view_g.el(2) = 1;
  view_g.el(3) = 0;
  
  //setcolor(140);   // za create_gif
  

  vector view_l(4);
  view_l = l_from_g * view_g;
 
  // form discrete arrays of eta and omega values to draw the "edges" correctly

  double eta_array[100];
  int eta_len = 0;
  

  double eta_up = atan2(a3, a2);
  double eta_down = -PI/2.0 + eta_up;
  double eta_over;

  for(eta = -PI/2.0 ; eta < eta_down; eta += PI/12.0)
    eta_array[eta_len++] = eta;
  eta_array[eta_len++] = eta_down;
  eta_over = eta;
  for(eta = eta_over; eta < eta_up; eta += PI/12.0)
    eta_array[eta_len++] = eta;
  eta_array[eta_len++] = eta_up;
  eta_over = eta;
  for(eta = eta_over; eta < PI/2.0; eta += PI/12.0)
    eta_array[eta_len++] = eta;
  eta_array[eta_len++] = PI/2.0;

  double omega_array[100];
  int omega_len = 0;

  for(omega = -PI; omega < +PI; omega += PI/12.0)
    omega_array[omega_len++] = omega;
  omega_array[omega_len++] = +PI;


  // draw const eta lines

  for(i = 0; i < eta_len; i++) {
 
    eta = eta_array[i];
    j = 0;
    omega = omega_array[j];
    old_r = g_from_l*(this->r(map_eta(eta), map_omega(omega)));

    for(j = 1; j < omega_len; j++) {

      omega = omega_array[j];
      new_r = g_from_l*((*this).r(map_eta(eta), map_omega(omega)));

      // draw the line if visible or no_hidding = TRUE

      if (no_hidding || ((view_l * normal(map_eta(eta), map_omega(omega))) >= -0.1))
      { old_r1 = transform(old_r);
        new_r1 = transform(new_r);
        //line((int) old_r1.el(0), ih() - (int) old_r1.el(1), (int) new_r1.el(0), ih() - (int) new_r1.el(1));
        }
        
        old_r = new_r;
    }

    //flush();
 
  }

  // draw const omega lines

  for(i = 0; i < omega_len; i++) {
 
    omega = omega_array[i];
    j = 0;
    eta = eta_array[j];
    old_r = g_from_l*(this->r(map_eta(eta), map_omega(omega)));

    for(j = 1; j < eta_len; j++) {

      eta = eta_array[j];
      new_r = g_from_l*((*this).r(map_eta(eta), map_omega(omega)));

      // draw the line if visible

      if (no_hidding || (view_l * normal(map_eta(eta), map_omega(omega))) >= -0.1)
      { old_r1 = transform(old_r);
        new_r1 = transform(new_r);
        //line((int) old_r1.el(0), ih() - (int) old_r1.el(1), (int) new_r1.el(0), ih() - (int) new_r1.el(1));        
        }
        old_r = new_r;
    }
  }

  //flush();  
}

void sq::print() 
{
  
  std::cout << "\nSuperellipsoid parameters:\n";
  std::cout << "a1 = " << a1 << " a2 = " << a2 << " a3 = " << a3 << '\n';
  std::cout << "e1 = " << e1 << " e2 = " << e2 << '\n';
  std::cout << "phi = " << phi << " theta = " << theta << " psi = " << psi << '\n';
  std::cout << "px = " << px << " py  = " << py << " pz = " << pz << '\n';  
  std::cout << g_from_l;
  std::cout << l_from_g; 

}  

void sq::fprint(FILE *f)
{ fprintf(f,"2 %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",a1,a2,a3,e1,e2,px,py,pz,theta,psi);
}

void sq::save(FILE *f)
{ fwrite(&a1, 1, sizeof(double), f);
  fwrite(&a2, 1, sizeof(double), f);
  fwrite(&a3, 1, sizeof(double), f);
  fwrite(&e1, 1, sizeof(double), f);
  fwrite(&e2, 1, sizeof(double), f);
  fwrite(&phi, 1, sizeof(double), f);
  fwrite(&theta, 1, sizeof(double), f);
  fwrite(&psi, 1, sizeof(double), f);
  fwrite(&px, 1, sizeof(double), f);
  fwrite(&py, 1, sizeof(double), f);
  fwrite(&pz, 1, sizeof(double), f);
  l_from_g.save(f);
  g_from_l.save(f);
  }
  
model* sq::improve(region& orig, region& added)
{ region r = orig | added;
  return(new sq(r));
  }
  
model* sq::improve(model *m, region& orig, region& added)
{ region r = orig | added;
  return(new sq((sq *)m,r));
  }

vector sq::transform(vector& v)
{ vector v1(4);
  v1.el(0) = iw() * (v.el(0) - minx) / dx;
  v1.el(1) = ih() * (v.el(1) - miny) / dy;
  return(v1);
  }  

void sq::parameters(char **name, double *value)
{ sprintf(name[0],"a1"); value[0] = a1;
  sprintf(name[1],"a2"); value[1] = a2;
  sprintf(name[2],"a3"); value[2] = a3;
  sprintf(name[3],"e1"); value[3] = e1;
  sprintf(name[4],"e2"); value[4] = e2;
  sprintf(name[5],"phi"); value[5] = phi;
  sprintf(name[6],"theta"); value[6] = theta;
  sprintf(name[7],"psi"); value[7] = psi;
  sprintf(name[8],"px"); value[8] = px;
  sprintf(name[9],"py"); value[9] = py;
  sprintf(name[10],"pz"); value[10] = pz;

  }

void sq::set_parameters(double *value)
{ a1 = value[0];
  a2 = value[1];
  a3 = value[2];
  e1 = value[3];
  e2 = value[4];
  phi = value[5];
  theta = value[6];
  psi = value[7];
  px = value[8];
  py = value[9];
  pz = value[10];

  // TODO: check if needed
  int sing;
  g_from_l.identity();
  l_from_g.identity();
  
  g_from_l.translate_g(px, py, pz);
  g_from_l.rotate_zr(phi);
  g_from_l.rotate_yr(theta);
  g_from_l.rotate_zr(psi);
  
  l_from_g = g_from_l.inverse(sing);
  print();
  }
