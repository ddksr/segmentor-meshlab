#include "segDrawer.h"
#include "segMesh.h"
#include "libsegmentor/common.h"
#include "libsegmentor/model.h"
#include "libsegmentor/sq.h"
#include "libsegmentor/asq.h"
#include "libsegmentor/tsq.h"

#include <common/meshmodel.h>
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

  qDebug() << "Models:" << models.size() << "Regions:" << regions.size();
}

void MeshlabDrawer::draw() {
  if (glewInit() != GLEW_OK) {
	qDebug() << "Problem :(";
	return;
  }
  qDebug() << "Models:" << models.size() << "Regions:" << regions.size();
  for (int i = 0; i < models.size(); ++i) {
	model* m = models.at(i);
	Color* c = modelColors.at(i);
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
	case CASQ:
	  draw_asq((asq*) m, c);
	  break;
	case CTSQ:
	  draw_tsq((tsq*) m, c);
	  break;
	case CBSQ:
	  draw_sq((sq*) m, c);
	  break;
	}
  }
  for (int i = 0; i < regions.size(); ++i) {
	region* r = regions.at(i);
	Color* c = regionColors.at(i);
	//draw_region(r, c);
  }
}

void MeshlabDrawer::draw_region(region* r, Color* c) {
  int i,j,n,n1,k,l,j1,j2,rob;
  struct point p,p1;
  static int cl = 0;

  glPushMatrix();
  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glEnable(GL_COLOR_MATERIAL);
  
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
			glVertex3f(p.x,p.y,p.z);
			p = r->theImage->pixel(j1);
			glVertex3f(p.x,p.y,p.z);
			p = r->theImage->pixel(j2);
			glVertex3f(p.x,p.y,p.z);   
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
				  glVertex3f(p.x,p.y,p.z);
				  glVertex3f(p1.x,p1.y,p1.z);
				  glEnd();
				}
			  } 
			}
		  }
		}
	  }
	}
  }
  glFlush();
  glPopMatrix();
  glPopAttrib();
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glMatrixMode(GL_MODELVIEW);
}

void MeshlabDrawer::draw_sq(sq* s, Color* c) {
  double eta,omega;
  vector vert(4);

  glPushMatrix();
  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glEnable(GL_COLOR_MATERIAL);
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  
  for (eta = -PI/2; eta <= PI/2; eta+=PI/12)         // constant eta lines
	{ glBegin(GL_POLYGON);
	  for (omega = -PI; omega <= PI; omega+=PI/12)
		{ vert = s->g_from_l*s->r(s->map_eta(eta),s->map_omega(omega));
		  glColor3f(c->r, c->g, c->b);
		  glVertex3f(vert.el(0),vert.el(1),vert.el(2));
		}
	  glEnd();
    } 
  for (omega = -PI/2; omega <= PI/2; omega+=PI/12)         // constant omega lines
	{ glBegin(GL_POLYGON);
	  
	  for (eta = -PI; eta <= PI; eta+=PI/12)
		{ vert = s->g_from_l*s->r(s->map_eta(eta),s->map_omega(omega));
		  glColor3f(c->r, c->g, c->b);
		  glVertex3f(vert.el(0),vert.el(1),vert.el(2));
		}
	  glEnd();
    } 

  glFlush();
  glPopMatrix();
  glPopAttrib();
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glMatrixMode(GL_MODELVIEW);
}

void MeshlabDrawer::draw_asq(asq* s, Color* c) {
  double eta,omega;
  vector vert(4);

  glPushMatrix();
  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glEnable(GL_COLOR_MATERIAL);
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  double step = s->kf > 0.1 ? PI/24 : PI/12;
  
  for (eta = -PI/2; eta <= PI/2; eta+=step)         // constant eta lines
	{ glBegin(GL_POLYGON);
	  for (omega = -PI; omega <= PI; omega+=step)
		{ vert = s->g_from_l*s->r(s->map_eta(eta),s->map_omega(omega));
		  glColor3f(c->r, c->g, c->b);
		  glVertex3f(vert.el(0),vert.el(1),vert.el(2));
		}
	  glEnd();
    } 
  for (omega = -PI/2; omega <= PI/2; omega+=step)         // constant omega lines
	{ glBegin(GL_POLYGON);
	  
	  for (eta = -PI; eta <= PI; eta+=step)
		{ vert = s->g_from_l*s->r(s->map_eta(eta),s->map_omega(omega));
		  glColor3f(c->r, c->g, c->b);
		  glVertex3f(vert.el(0),vert.el(1),vert.el(2));
		}
	  glEnd();
    } 

  glFlush();
  glPopMatrix();
  glPopAttrib();
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glMatrixMode(GL_MODELVIEW);
}

void MeshlabDrawer::draw_tsq(tsq* s, Color* c) {
  qDebug() << "draw tsq";
  double eta,omega;
  vector vert(4);

  glPushMatrix();
  glPushAttrib(GL_ALL_ATTRIB_BITS);
  glEnable(GL_COLOR_MATERIAL);
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  
  for (eta = -PI/2; eta <= PI/2; eta+=PI/12)         // constant eta lines
	{ glBegin(GL_POLYGON);
	  for (omega = -PI; omega <= PI; omega+=PI/12)
		{ vert = s->g_from_l*s->r(s->map_eta(eta),s->map_omega(omega));
		  glColor3f(c->r, c->g, c->b);
		  glVertex3f(vert.el(0),vert.el(1),vert.el(2));
		}
	  glEnd();
    } 
  for (omega = -PI/2; omega <= PI/2; omega+=PI/12)         // constant omega lines
	{ glBegin(GL_POLYGON);
	  
	  for (eta = -PI; eta <= PI; eta+=PI/12)
		{ vert = s->g_from_l*s->r(s->map_eta(eta),s->map_omega(omega));
		  glColor3f(c->r, c->g, c->b);
		  glVertex3f(vert.el(0),vert.el(1),vert.el(2));
		}
	  glEnd();
    } 

  glFlush();
  glPopMatrix();
  glPopAttrib();
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glMatrixMode(GL_MODELVIEW);
}

void MeshlabDrawer::draw_plane(plane* m, Color* c) {
  return;
  MeshModel* model = ((segMesh*)m->theImage)->mesh;
  CMeshO::FaceIterator fi;

  glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
  
  for(fi=model->cm.face.begin(); fi!=model->cm.face.end(); ++fi) {
	vcg::Point3f ma = fi->cV(0)->P(), mb = fi->cV(1)->P(), mc = fi->cV(2)->P();
	struct point a, b, c;
	
	a.x = ma[0];
	a.y = ma[1];
	a.z = ma[2];

	b.x = mb[0];
	b.y = mb[1];
	b.z = mb[2];
	
	c.x = mc[0];
	c.y = mc[1];
	c.z = mc[2];

	double dist = abs(m->abs_signed_distance(a));
	dist += abs(m->abs_signed_distance(b));
	dist += abs(m->abs_signed_distance(c));

	if (dist < 0.1) {
	  // (*fi).C()[0] = 255;
	  // (*fi).C()[0] = 0;
	  // (*fi).C()[0] = 0;
	  // (*fi).C()[0] = 255;
	}
	
  }  
}
void MeshlabDrawer::draw_surface2(surface2* m, Color* c) {
  // TODO
}

void MeshlabDrawer::draw_cone(cone* m, Color* c) {
  glPushMatrix();
  glPushAttrib(GL_ALL_ATTRIB_BITS);
  
  vector3D axis,n, orig, n2, v1, v2, point, m1, m2, axis2, min_axis, axis3, axis4;
  vector3D axispoint;

  int rot = (int)fabs(1/m->a[6]), i, j;
  double dz;
  double kot, rad, pi = 3.1415926535, rad1;

  axis.el(0) = cos(m->a[4]) * sin(m->a[5]);
  axis.el(1) = sin(m->a[4]) * sin(m->a[5]);
  axis.el(2) = cos(m->a[5]);

  min_axis = axis;
  if (m->a[6] > 0.0) min_axis.multiply_col(0,-1);

  n.el(0) = cos(m->a[2]) * sin(m->a[3]);
  n.el(1) = sin(m->a[2]) * sin(m->a[3]);
  n.el(2) = cos(m->a[3]);

  kot = acos(n.inner(min_axis));
  axis3 = axis;
  axis3.multiply_col(0,(1/fabs(m->a[6]))/cos(kot));

  for (i = 0; i < 3; i++) orig.el(i) = m->a[11+i];

  n2 = n;
  n2.multiply_col(0,m->a[1]+1/m->a[6]);

  v1 = axis.outer(n);
  v1 = v1.normalize();
  v2 = axis.outer(v1);
  
  if (rot > 100) rot = 100;
  if (rot < 30) rot = 30;

  glColor3f(0.0,1.0,0.0);
  glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);  
  
  dz = (m->zmax-m->zmin)/10;

  for (j = 0; j <= 10; j++)
  { rad = (m->zmin + j * dz) / tan(kot);
    axis2 = axis;
    axis2.multiply_col(0,m->zmin+j*dz);
     
    glBegin(GL_POLYGON);
    axispoint = orig + n2 + axis3 + axis2;
    for (i = 0; i < rot; i++)
    { m1 = v1;
      m2 = v2;
      
      m1.multiply_col(0,rad * cos(2 * i * pi / rot));
      m2.multiply_col(0,rad * sin(2 * i * pi / rot));
      point = axispoint + m1 + m2;
      glVertex3f(point.el(0),point.el(1),point.el(2));
      }
    glEnd();
    } 

  axis4 = orig + n2 + axis3;

  rad = m->zmin / tan(kot);
  rad1 = m->zmax / tan(kot);

  axis2 = axis;
  axis2.multiply_col(0,m->zmin);
  axis4 = axis;
  axis4.multiply_col(0,m->zmax);

  for (i = 0; i < rot; i++)
  { glBegin(GL_LINES);
    m1 = v1;
    m2 = v2;
    m1.multiply_col(0,rad1 * cos(2 * i * pi / rot));
    m2.multiply_col(0,rad1 * sin(2 * i * pi / rot));    
    point = orig + n2 + axis3 + axis4 + m1 + m2;
    glVertex3f(point.el(0),point.el(1),point.el(2));
    m1 = v1;
    m2 = v2;
    m1.multiply_col(0,rad * cos(2 * i * pi / rot));
    m2.multiply_col(0,rad * sin(2 * i * pi / rot));
    point = orig + n2 + axis3 + axis2 + m1 + m2;
    glVertex3f(point.el(0),point.el(1),point.el(2));

    glEnd();
    }

  glFlush();
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);


}

void MeshlabDrawer::draw_torus(torus* m, Color* col) {
  glPushMatrix();
  glPushAttrib(GL_ALL_ATTRIB_BITS);
  
  vector3D axis,n,r1,r2,c,n2,a2,m1,m2,r,rm,m3,m4,pt;

  double rot1, rot2;
 
  double pi = 3.1415926535;
  double kot,h,vr;
  int i,j;

  axis.el(0) = cos(m->a[4])*sin(m->a[5]);
  axis.el(1) = sin(m->a[4])*sin(m->a[5]);
  axis.el(2) = cos(m->a[5]);

  n.el(0) = cos(m->a[2])*sin(m->a[3]);
  n.el(1) = sin(m->a[2])*sin(m->a[3]);
  n.el(2) = cos(m->a[3]);

  std::cout << "Axis:" << axis << "N:" << n;

  kot = acos(n.inner(axis));

  r1 = axis.outer(n);
  r1 = r1.normalize();
  r2 = axis.outer(r1);  
  r2 = r2.normalize();

  n2 = n;
  n2.multiply_col(0,m->a[1]+1/m->a[7]);

  a2 = axis;
  h = fabs(1/m->a[7]-1/m->a[6])*cos(kot);
  a2.multiply_col(0,h);
  vr = fabs(1/m->a[7]-1/m->a[6])*sin(kot);  

  c = n2 + a2;

  c.el(0) += m->a[11];
  c.el(1) += m->a[12];
  c.el(2) += m->a[13];

  rot2 = 1/m->a[6];
  rot1 = vr;

  if (rot1 < 30) rot1 = 30;
  if (rot1 > 100) rot1 = 100;
  
  if (rot2 < 20) rot2 = 20;
  if (rot2 > 50) rot2 = 50;

  rm = r1;
  rm.multiply_col(0,vr);

  glColor3f(0.0,0.0,1.0);
  glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);  

  for (i = 0; i < rot1; i++)
  { glBegin(GL_POLYGON);

    m1 = r1;
    m2 = r2;
    m1.multiply_col(0,cos(i*2*pi/rot1));
    m2.multiply_col(0,sin(i*2*pi/rot1));
    r = m1 + m2;
    rm = r;
    rm.multiply_col(0,vr);

    for (j = 0; j < rot2; j++)
    { m3 = axis;
      m4 = r;
      m3.multiply_col(0,cos(j*2*pi/rot2)/m->a[6]);
      m4.multiply_col(0,sin(j*2*pi/rot2)/m->a[6]);
      pt = c + rm + m3 + m4;
      glVertex3f(pt.el(0),pt.el(1),pt.el(2));
      }
    glEnd();
    }

  for (j = 0; j < rot2; j++)
  { glBegin(GL_POLYGON);
    m3 = axis;
    m3.multiply_col(0,cos(j*2*pi/rot2)/m->a[6]);
    
    for (i = 0; i < rot1; i++)
    { m1 = r1;
      m2 = r2;
      m1.multiply_col(0,cos(i*2*pi/rot1));
      m2.multiply_col(0,sin(i*2*pi/rot1));
      r = m1 + m2;
      rm = r;
      rm.multiply_col(0,vr);
      m4 = r.normalize();
      m4.multiply_col(0,sin(j*2*pi/rot2)/m->a[6]);
      pt = c + rm + m3 + m4;  
      glVertex3f(pt.el(0),pt.el(1),pt.el(2));
      } 

    glEnd();
    } 

  glFlush();
  glPopMatrix();
  glPopAttrib();
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glMatrixMode(GL_MODELVIEW);
}

void MeshlabDrawer::draw_sphere(sphere* m, Color* c) {
  glPushMatrix();
  glPushAttrib(GL_ALL_ATTRIB_BITS);

  double ksi, phi;
  double dk = PI / 12;
  double va,vb,vc;

  glColor3f(0.0,0.0,0.0);

  for (phi = -PI / 2; phi < PI / 2; phi+=dk)     // const phi circles
  { glBegin(GL_LINE_LOOP);
    for (ksi = -PI; ksi < PI; ksi+=dk)
    { va = m->a + m->radius * cos(phi) * cos(ksi);
      vb = m->b + m->radius * cos(phi) * sin(ksi);
      vc = m->c + m->radius * sin(phi);
      // dist = sqrt((va-a)*(va-a)+(vb-b)*(vb-b)+(vc-c)*(vc-c)) - m->radius;
      glVertex3f(va,vb,vc);
      }
    glEnd();
    }

  for (ksi = -PI / 2; ksi < PI / 2; ksi+=dk)            // const ksi circles
  { glBegin(GL_LINE_LOOP);
    for (phi = -PI; phi < PI; phi+=dk)
    { va = m->a + m->radius * cos(phi) * cos(ksi);
      vb = m->b + m->radius * cos(phi) * sin(ksi);
      vc = m->c + m->radius * sin(phi);
      // dist = sqrt((va-a)*(va-a)+(vb-b)*(vb-b)+(vc-c)*(vc-c)) - m->radius;
      glVertex3f(va,vb,vc);
      }
    glEnd();
    }

  
  glFlush();
  glPopMatrix();
  glPopAttrib();
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glMatrixMode(GL_MODELVIEW);
}

void MeshlabDrawer::draw_cylinder(cylinder* m, Color* c) {
  glPushMatrix();
  glPushAttrib(GL_ALL_ATTRIB_BITS);

  vector3D axis,n,n1,n2,m1,m2,vr,axis2;
  int i,j,rot = (int)(1/fabs(m->a[5]));
  double d = m->a[1] + 1 / m->a[5], z;
  double dz = (m->zmax - m->zmin)/10;
  if (rot > 100) rot = 100;
  if (rot < 30) rot = 30;

  glColor3f(1.0,0.0,0.0);
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  axis.el(0) = cos(m->a[2])*cos(m->a[3])*cos(m->a[4]) - sin(m->a[2])*sin(m->a[4]);
  axis.el(1) = sin(m->a[2])*cos(m->a[3])*cos(m->a[4]) + cos(m->a[2])*sin(m->a[4]);
  axis.el(2) = -sin(m->a[3])*cos(m->a[4]);
    
  n.el(0) = cos(m->a[2])*sin(m->a[3]);
  n.el(1) = sin(m->a[2])*sin(m->a[3]);
  n.el(2) = cos(m->a[3]);

  n1 = n;

  n.multiply_col(0,d);

  n2 = axis.outer(n1);

  for (j = 0; j <= 10; j++)
  { z = m->zmin + j * dz;
    glBegin(GL_POLYGON);
    if (j == 10) glColor3f(0.0,0.0,1.0);
    for (i = 0; i < rot; i++)
    { m1 = n1;
      m2 = n2;
      axis2 = axis;
      m1.multiply_col(0,cos(2 * i * PI / rot) / m->a[5]);
      m2.multiply_col(0,sin(2 * i * PI / rot) / m->a[5]);
      axis2.multiply_col(0,z);
      vr = n + axis2 + m1 + m2;
      glVertex3f(vr.el(0)+m->a[11],vr.el(1)+m->a[12],vr.el(2)+m->a[13]);
      }
    glEnd();
    }

  glColor3f(1.0,0.0,0.0);
  if (rot % 2) rot++;
  for (i = 0; i < rot; i++)
  { glBegin(GL_QUADS);
    m1 = n1;
    m2 = n2;
    axis2 = axis;
    m1.multiply_col(0,cos(2 * i * PI / rot) / m->a[5]);
    m2.multiply_col(0,sin(2 * i * PI / rot) / m->a[5]);
    axis2.multiply_col(0,m->zmin);
    vr = n + axis2 + m1 + m2;
    glVertex3f(vr.el(0)+m->a[11],vr.el(1)+m->a[12],vr.el(2)+m->a[13]);
    axis2 = axis;
    axis2.multiply_col(0,m->zmax);
    vr = n + axis2 + m1 + m2;
    glVertex3f(vr.el(0)+m->a[11],vr.el(1)+m->a[12],vr.el(2)+m->a[13]);

    m1 = n1;
    m2 = n2;
    axis2 = axis;
    m1.multiply_col(0,cos(PI + 2 * i * PI / rot) / m->a[5]);
    m2.multiply_col(0,sin(PI + 2 * i * PI / rot) / m->a[5]);
    axis2.multiply_col(0,m->zmax);
    vr = n + axis2 + m1 + m2;
    glVertex3f(vr.el(0)+m->a[11],vr.el(1)+m->a[12],vr.el(2)+m->a[13]);
    axis2 = axis;
    axis2.multiply_col(0,m->zmin);
    vr = n + axis2 + m1 + m2;
    glVertex3f(vr.el(0)+m->a[11],vr.el(1)+m->a[12],vr.el(2)+m->a[13]);
    glEnd();
    }
  

  glFlush();
  glPopMatrix();
  glPopAttrib();
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glMatrixMode(GL_MODELVIEW);
}
