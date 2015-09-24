#include "segDrawer.h"
#include "libsegmentor/common.h"
#include "libsegmentor/model.h"
#include "libsegmentor/sq.h"

#include <common/interfaces.h>

#include <GL/glew.h>

// using namespace vcg;

MeshlabDrawer::MeshlabDrawer(image* im) {
  img = im;
}

void MeshlabDrawer::prepare(model* m) {
  prepare(m, 1.0, 0.0, 0.0);
}

void MeshlabDrawer::prepare(model* m, float r, float g, float b) {
  models << m;
  modelColors << new Color{r, g, b};
}


void MeshlabDrawer::prepare(region* r) {
  prepare(r, 1.0, 0.0, 0.0);
}

void MeshlabDrawer::prepare(region* reg, float r, float g, float b) {
  regions << reg;
  regionColors << new Color{r, g, b};
}

void MeshlabDrawer::clear() {
  models.clear();
  modelColors.clear();
  regions.clear();
  regionColors.clear();
}

void MeshlabDrawer::draw() {
  if (glewInit() != GLEW_OK) {
	qDebug() << "Problem :(";
	return;
  }
  for (int i = 0; i < models.size(); ++i) {
	model* m = models.at(i);
	Color* c = modelColors.at(i);
	qDebug() << "Drawing model" << m->what_model();
	switch(m->what_model()) {
	case CPLANE:
	  draw_plane((plane*) m, c);
	  break;
	case CSQ:
	  draw_sq((sq*) m, c);
	  break;
	case CSURFACE2:
	  draw_surface2((surface2*) m, c);
	  break;
	case CSPHERE:
	  draw_sphere((sphere*) m, c);
	  break;
	case CCYLINDER:
	  draw_cylinder((cylinder*) m, c);
	  break;
	case CCONE:
	  draw_cone((cone*) m, c);
	  break;
	case CTORUS:
	  draw_torus((torus*) m, c);
	  break;
	}
  }
  for (int i = 0; i < regions.size(); ++i) {
	qDebug() << "Drawing region";
	region* r = regions.at(i);
	Color* c = regionColors.at(i);
  }
}

void MeshlabDrawer::draw_region(region* r, Color* c) {
  int i,j,n,n1,k,l,j1,j2,rob;
  struct point p,p1;
  static int cl = 0;
  
  cl = (cl + 1) % 12;
 
  for (j = r->miny; j <= r->maxy; j++) {
    for (i = r->minx; i <= r->maxx; i++) {
      if (r->included(i,j)) {
		n = r->theImage->neighbours(i,j)>>1;
		rob = 0;
		for (k = 0; k < n; k++) {
		  j1 = r->theImage->neighbour(i,0,k<<1);
		  j2 = r->theImage->neighbour(i,0,1+(k<<1));
		  if (r->included(j1) && r->included(j2)) {
			if (!(cl%3)) glColor3f(0.5,0.6+(cl/3)*0.1,0.6+(cl/3)*0.1);
			else if (cl%3 == 1) glColor3f(0.6+(cl/3)*0.1,0.5,0.6+(cl/3)*0.1);
			else glColor3f(0.6+(cl/3)*0.1,0.6+(cl/3)*0.1,0.5);
			  
			glBegin(GL_TRIANGLES);
			p = r->theImage->pixel(i,j);
			glVertex4f(p.x,p.y,p.z,r->theImage->norm);
			//if (r->theImage->nGL) glNormal3f(p.nx,p.ny,p.nz);
			p = r->theImage->pixel(j1);
			glVertex4f(p.x,p.y,p.z,r->theImage->norm);
			//if (r->theImage->nGL) glNormal3f(p.nx,p.ny,p.nz);
			p = r->theImage->pixel(j2);
			glVertex4f(p.x,p.y,p.z,r->theImage->norm);   
			//if (r->theImage->nGL) glNormal3f(p.nx,p.ny,p.nz);         
			glEnd();
		  } else { rob = 1; }
		  if (rob) {
			n *= 2;
			p = r->theImage->pixel(i,j);
			for (k = 0; k < n; k++){
			  j1 = r->theImage->neighbour(i,0,k);
			  if (r->included(j1)) {
				p1 = r->theImage->pixel(j1,0); 
				n1 = r->theImage->neighbours(j1,0);              
				rob = 0;
				for (l = 0; l < n1; l++) {
				  j2 = r->theImage->neighbour(j1,0,l);
				  if (!r->included(j2)) rob = 1;
				}
				if (rob) {
				  glBegin(GL_LINES);
				  glColor3f(0.0,0.0,0.0);
				  glVertex4f(p.x,p.y,p.z,r->theImage->norm);
				  glVertex4f(p1.x,p1.y,p1.z,r->theImage->norm);
				  glEnd();
				}
			  } 
			}
		  }
		}
	  }
	}
  }
}

void MeshlabDrawer::draw_sq(sq* s, Color* c) {
  double eta,omega;
  vector vert(4);
  
  glPushMatrix();
  glColor3f(c->r, c->g, c->b);
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    
  for (eta = -PI/2; eta <= PI/2; eta+=PI/12)         // constant eta lines
	{ glBegin(GL_POLYGON);
	  for (omega = -PI; omega <= PI; omega+=PI/12)
		{ vert = s->g_from_l*s->r(s->map_eta(eta),s->map_omega(omega));
		  glVertex4f(vert.el(0),vert.el(1),vert.el(2),img->norm);
		}
	  glEnd();
    } 
  for (omega = -PI/2; omega <= PI/2; omega+=PI/12)         // constant omega lines
	{ glBegin(GL_POLYGON);
	  for (eta = -PI; eta <= PI; eta+=PI/12)
		{ vert = s->g_from_l*s->r(s->map_eta(eta),s->map_omega(omega));
		  glVertex4f(vert.el(0),vert.el(1),vert.el(2),img->norm);
		}
	  glEnd();
    } 

  glFlush();
  glPopMatrix();
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glMatrixMode(GL_MODELVIEW);
}

void MeshlabDrawer::draw_plane(plane* m, Color* c) {}
void MeshlabDrawer::draw_surface2(surface2* m, Color* c) {}
void MeshlabDrawer::draw_cone(cone* m, Color* c) {}
void MeshlabDrawer::draw_torus(torus* m, Color* c) {}
void MeshlabDrawer::draw_sphere(sphere* m, Color* c) {}
void MeshlabDrawer::draw_cylinder(cylinder* m, Color* c) {}
