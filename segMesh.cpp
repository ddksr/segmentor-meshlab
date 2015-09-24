#include <Qt>

#include <common/interfaces.h>

#include "libsegmentor/image.h"
#include "libsegmentor/plane.h"
#include "libsegmentor/region.h"
#include "segMesh.h"

using namespace vcg;

segMesh::segMesh(MeshModel *model) : image() {
  CMeshO::VertexIterator vi;
  CMeshO::FaceIterator fi;
  int i = 0;
  double x1, y1, z1, dx, dy, dz;
  for(vi=model->cm.vert.begin();vi!=model->cm.vert.end();++vi) {
	Point3f p = vi->P();
	x1 = (double)p[0];
	y1 = (double)p[1];
	z1 = (double)p[2];

	if (!i){
	  minx = maxx = x1;
      miny = maxy = y1;
      minz = maxz = z1;
	} else {
	  if (x1 < minx) minx = x1;
      if (x1 > maxx) maxx = x1;
      if (y1 < miny) miny = y1;
      if (y1 > maxy) maxy = y1;
      if (z1 < minz) minz = z1;
      if (z1 > maxz) maxz = z1;
	}
	i++;
  }
  dx = maxx - minx;
  dy = maxy - miny;
  dz = maxz - minz;
  norm = max(dx,max(dy,dz));
  qDebug("Points: %d", i);
  pix = new struct point[number = w = i];
  h = 1;
  i = 0;
  x = new int[number];
  y = new int[number];
  for(vi=model->cm.vert.begin();vi!=model->cm.vert.end();++vi) {
	Point3f p = vi->P();
	pix[i].x = (double)p[0];
	pix[i].y = (double)p[1];
	pix[i].z = (double)p[2];
	pix[i].n = 0;
	pix[i].valid = 1;
	pix[i].boundary = NO_BOUND;
	
    x[i] = (int)(256 * (pix[i].x - minx) / dx);
    y[i] = (int)(256 * (pix[i].y - miny) / dy);

	i++;
  }

  i = 0;

  const CMeshO::VertexType * v0 = &(model->cm.vert[0]);
  for(fi=model->cm.face.begin();fi!=model->cm.face.end();++fi) {
	uint a = int(fi->cV(0) - v0);
	uint b = int(fi->cV(1) - v0);
	uint c = int(fi->cV(2) - v0);

	add_neighbour(a,0,b);
	add_neighbour(a,0,c);
	add_neighbour(b,0,a);
	add_neighbour(b,0,c);
	add_neighbour(c,0,a);
	add_neighbour(c,0,b);

	i++;
  }
  qDebug("Triangles: %d", i);  
}

segMesh::segMesh(int w1, int h1) : image()
{
  unsigned int i;
  w = w1;
  h = h1;
  pix = new struct point[number = w * h];
  for (i = 0; i < number; i++) 
  { pix[i].n = 0;
    pix[i].valid = 1;
    }
  x = NULL;
  y = NULL;
}  


segMesh::~segMesh() {}

void segMesh::bounding_box(double &x1, double &y1, double &z1, double &x2, double &y2, double &z2) {
  x1 = minx;
  y1 = miny;
  z1 = minz;
  x2 = maxx;
  y2 = maxy;
  z2 = maxz;
}

void segMesh::extremes() {
  for (int i = 0; i < number; i++) {
	if (!i || pix[i].x < minx) minx = pix[i].x;
    if (!i || pix[i].x > maxx) maxx = pix[i].x;
    if (!i || pix[i].y < miny) miny = pix[i].y;
    if (!i || pix[i].y > maxy) maxy = pix[i].y;
    if (!i || pix[i].z < minz) minz = pix[i].z;
    if (!i || pix[i].z > maxz) maxz = pix[i].z;
  }
}

void segMesh::mark_boundary()
{ int i,j,k,num;
  int type[4];
  int n2 = number;

  for (i = 0; i < n2; i++)
  { if (!pix[i].n) pix[i].boundary = OUTLIER;
    if (pix[i].boundary == NO_BOUND)
    { type[1] = type[2] = type[3] = 0;
      for (j = 0; j < pix[i].n; j++)
      { num = 0; 
        for (k = 0; k < pix[i].n; k++)
          if (pix[i].neigh[k] == pix[i].neigh[j]) num++;
        if (num > 2) type[3] = 1;
          else type[num] = 1;
        }
      if (type[3]) pix[i].boundary = IRREG;
        else if (type[1]) pix[i].boundary = BOUND;
      }      
    }
  }

image* segMesh::calcNormals(){
  int i,w1,h1,n;
  double x,y,z;
  region r(this);
  plane *p;
  segMesh *c_normals;

  c_normals = new segMesh(w1 = this->width(),h1 = this->height());

  n = w1 * h1;
  for (i = 0; i < n; i++)
  { r.reset_all();
    r.set_point(i);
    r |= r.neighbourhood();
    r |= r.neighbourhood();
    if (r.point_count() >= 3) 
    { p = new plane(r);
      p -> get_normal_vector(x,y,z);
      c_normals -> set_pixel(i,0,x,y,z); 
      this->set_normal(i,x,y,z);
      delete p;
      }  
    }
  return (image*)c_normals;
}  



//segMesh& segMesh::operator=(const segMesh &im) {}
