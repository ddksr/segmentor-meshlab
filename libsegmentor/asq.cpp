// sq.C

#include "asq.h"
#include "state.h"
#include "mathutil.h"

#include <math.h>
#include <iostream>

extern "C" {
#include "sq_recover.h"
 double sq_mindist(double a1, double a2, double a3, double e1, double e2, double x0, double y0, double z0);
}

int recover_params(asq* model, vect* list, int no) {
  double _bk= 0.0000010, _ba = 0, _kx = 0, _ky = 0;
  if (State::lookup->on && State::lookup->isSet) {
	// Reset SQ params to those from lookup
	model->a1 = State::lookup->a;
	model->a2 = State::lookup->b;
	model->a3 = State::lookup->c;
	model->e1 = State::lookup->e1;
	model->e2 = State::lookup->e2;
	
	return recover_search(list, no,
						  &model->a1, &model->a2, &model->a3, &model->e1, &model->e2,
						  &model->px, &model->py, &model->pz, &model->phi, &model->theta, &model->psi,
						  &_kx, &_ky, &_bk, &_ba, &model->kxy, &model->kf, RECOVER_ASQ);
  }
  // return recover(list, no,
  // 				  &model->a1, &model->a2, &model->a3, &model->e1, &model->e2,
  // 				  &model->px, &model->py, &model->pz, &model->phi, &model->theta, &model->psi);
  return recover2(list, no,
  				  &model->a1, &model->a2, &model->a3, &model->e1, &model->e2,
  				  &model->px, &model->py, &model->pz, &model->phi, &model->theta, &model->psi,
  				  &_kx, &_ky, &_bk, &_ba, &model->kxy, &model->kf, RECOVER_ASQ);
}

asq::asq(region &r) {

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

  kxy = 0;
  kf = 0;
  
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

asq::asq(asq *m, region& r) {
  #define MAX_LEN 25000

  struct vect list[MAX_LEN];
  int x, y, min_x, min_y, max_x, max_y, no;
  struct point p;
  int sing;
  
  a1 = m->a1; a2 = m->a2; a3 = m->a3; e1 = m->e1; e2 = m->e2;
  phi = m->phi; theta = m->theta; psi = m->psi; px = m->px; py = m->py; 
  pz = m->pz;

  kxy = m->kxy;
  kf = m->kf;
  
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


double asq::map_eta(double eta) {
  // TODO: maybe?
  double x = a2 * spow(cos(eta), 1/e1);
  double y = a3 * spow(sin(eta), 1/e1);
  
  return (atan2(y, x));
}

double asq::map_omega(double omega) {
  // TODO: maybe?
  double x = a2 * spow(cos(omega), 1/e2);
  double y = a1 * spow(sin(omega), 1/e2);

  return (atan2(y, x));
}

vector asq::normal(double eta, double omega) const {
  // TODO: maybe?
  vector result(4);

  result.el(0) = spow(cos(eta), 2.0 - e1) * spow(cos(omega), 2.0 - e2) / a1;
  result.el(1) = spow(cos(eta), 2.0 - e1) * spow(sin(omega), 2.0 - e2) / a2;
  result.el(2) = spow(sin(eta), 2.0 - e1) / a3;
  result.el(3) = 1.0;
  
  result = (1.0 / result.norm()) * (matrix)result;


  return result;
}

vector asq::r(double eta, double omega) const {
  
  vector result(4);

  double z = a3 * spow(sin(eta), e1),
	fx = kxy * z / a3 + 1,
	fy = kxy * z / a3 + 1;

  result.el(0) = a1 * spow(cos(eta), e1) * spow(cos(omega), e2) * fx;
  result.el(1) = a2 * spow(cos(eta), e1) * spow(sin(omega), e2) * fy;
  result.el(2) = z;
  result.el(3) = 1.0;
  
  return result;
}

double asq::abs_signed_distance(struct point& p)
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
  
double asq::f(double x, double y, double z) const {
  

  // thid double pow is to prevent negative arguments to pow when the exponent is negative
  // do not try to optimize

  double fx = kxy * z / a3 + 1.0,
	fy = kxy * z / a3 + 1.0; 

  double a = pow(pow((x / (a1 * fx)), 2), 1/e2);
  double b = pow(pow((y / (a2 * fy)), 2), 1/e2);   
  double c = pow(pow((z / a3), 2), 1/e1);
  
  return pow(pow(a + b, e2/e1) + c, e1);
}

void asq::print() 
{
  
  std::cout << "\nSuperellipsoid parameters:\n";
  std::cout << "a1 = " << a1 << " a2 = " << a2 << " a3 = " << a3 << '\n';
  std::cout << "e1 = " << e1 << " e2 = " << e2 << '\n';
  std::cout << "phi = " << phi << " theta = " << theta << " psi = " << psi << '\n';
  std::cout << "px = " << px << " py  = " << py << " pz = " << pz << '\n';
  std::cout << "kxy = " << kxy << " kf  = " << kf << '\n';  
  std::cout << g_from_l;
  std::cout << l_from_g; 

}  

model* asq::improve(region& orig, region& added)
{ region r = orig | added;
  return(new asq(r));
  }
  
model* asq::improve(model *m, region& orig, region& added)
{ region r = orig | added;
  return(new asq((asq *)m,r));
  }

vector asq::transform(vector& v)
{ vector v1(4);
  v1.el(0) = iw() * (v.el(0) - minx) / dx;
  v1.el(1) = ih() * (v.el(1) - miny) / dy;
  return(v1);
  }  

void asq::parameters(char **name, double *value)
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
  sprintf(name[11],"kxy"); value[11] = kxy;
  sprintf(name[12],"kf"); value[12] = kf;
  }

void asq::set_parameters(double *value)
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
  kxy = value[11];
  kf = value[12];

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
