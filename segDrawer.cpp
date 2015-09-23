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
  models << m;
}

void MeshlabDrawer::clear() {
  models.clear();
}

void MeshlabDrawer::draw() {
  for (int i = 0; i < models.size(); ++i) {
	MeshlabDrawer::draw_sq((sq*)models.at(i));
  }
}

void MeshlabDrawer::draw_sq(sq* s) {
  double eta,omega;
  vector vert(4);

  if (glewInit() != GLEW_OK) {
	qDebug() << "Problem :(";
  }
  glPushMatrix();
  glColor3f(1.0,0.0,0.0);
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

void MeshlabDrawer::draw_test() {

}
