#ifndef SEGMESH_H
#define SEGMESH_H

#include <common/interfaces.h>
#include <meshlabplugins/edit_pickpoints/pickedPoints.h>

#include "libsegmentor/image.h"
#include "libsegmentor/region.h"

class segMesh : public image
{
 public:
  MeshModel* mesh;
  
  segMesh(MeshModel*);
  segMesh(int, int);
  ~segMesh();
  
  int triang_count(int*);
  int triang_count();
  void bounding_box(double &x1, double &y1, double &z1, double &x2, double &y2, double &z2);
  int valid_point(struct point& p) { return(p.valid); }
  void extremes();
  void mark_boundary();
  
  image* calcNormals();

  void setSelectedPoints();

  
  
  //segMesh& operator=(const segMesh &im);
  
};

#endif
