#include "segmentation.h"
#include "planar_d.h"
#include "sq_d.h"
#include "surface2_d.h"
#include "sphere_d.h"
#include "cylinder_d.h"
#include "cone_d.h"
#include "torus_d.h"
#include "asq_d.h"
#include "tsq_d.h"
#include "bsq_d.h"
#include "image.h"
#include "transformation.h"

#include <iostream>

// TMP TMP TMP
#include <QDebug>

Drawer* segmentation::drawer = NULL;
ProgressIndicator* segmentation::progress = NULL;
Messaging* segmentation::message = NULL;  

std::ostream& operator<<(std::ostream&s, struct point& p)
{ int i;
  s << '(' << p.x << ',' << p.y << ',' << p.z << ") - neighbours:";
  for (i = 0; i < p.n; i++) s << ' ' << p.neigh[i];
  return(s);
  }

void segmentation::throw_away()
{ int i,j;
  for (j = i = 0; i < n; i++)
  { handle[j] = handle[i];
    if (d[handle[i]] != NULL) j++;
    }
  n = j;  
  }  

segmentation::segmentation(int nsize, int dnsize)
{ int i;
  
  if (d == NULL) 
  { d = new description*[dsize = dnsize];
    dn = 0;
    for (i = 0; i < dsize; i++) d[i] = NULL;
    
    /*   int num = im -> points();
    double minx,miny,maxx,maxy;
    
    p = im -> pixel(0);
    minx = maxx = p.x; 
    miny = maxy = p.y;
    for (i = 1; i < num; i++)
    { p = im -> pixel(i);
      if (p.x < minx) minx = p.x;
      if (p.x > maxx) maxx = p.x;
      if (p.y < miny) miny = p.y;
      if (p.y > maxy) maxy = p.y;
      }
    model::minx = minx;
    model::miny = miny;
    model::dx = maxx - minx;
    model::dy = maxy - miny; */ 
    }
  segmentationImage = NULL;
  normals = NULL;
  usingList++;
  n = 0;
  handle = new int[size = nsize];
  }

segmentation::segmentation(const segmentation& l)
{ int i;
  usingList++;
  handle = new int[size = l.size];
  segmentationImage = l.segmentationImage;
  normals = l.normals;
  n = l.n;
  k2 = l.k2; k3 = l.k3;
  for (i = 0; i < n; i++) handle[i] = l.handle[i];
  }  
  
segmentation& segmentation::operator=(const segmentation& l)
{ int i;
  if (this != &l)
  { delete [] handle;
    handle = new int[size = l.size];
    segmentationImage = l.segmentationImage;
    normals = l.normals;
    n = l.n;
    k2 = l.k2; k3 = l.k3;
    for (i = 0; i < n; i++) handle[i] = l.handle[i];
    }
  return(*this);  
  }  

segmentation& segmentation::operator+=(const segmentation& l)
{ int i;
  int new_size, *temp;
  if (this != &l && segmentationImage == l.segmentationImage)
  { temp = new int[new_size = size + l.size];
    for (i = 0; i < n; i++) temp[i] = handle[i];
    for (i = 0; i < l.n; i++) temp[i+n] = l.handle[i];
    delete [] handle;
    handle = temp;
    size += l.size;
    n += l.n;
    }
  return(*this);
  }
    
segmentation segmentation::operator+(const segmentation l)
{ segmentation lnew(*this);
  lnew += l;
  if (segmentationImage != l.segmentationImage) std::cout << "\nNot the same segmentationImage in operator +!\n";
  return(lnew);
  }

void segmentation::place_seeds(int seed_size, MODELTYPE type)
{ region r(segmentationImage), *rn, nb(segmentationImage);
  int i,j;
  
  for (j = 0; j < segmentationImage->height(); j++)
    for (i = 0; i < segmentationImage->width(); i++)
      if (segmentationImage->valid_point(r.get_point(i,j))) r.set_point(i,j);
  
  for (j = r.miny; j <= r.maxy; j++)
    for (i = r.minx; i <= r.maxx; i++) {
	  if (!r.is_selected(i, j)) continue;
      if (r.included(i,j)) {
		rn = new region(segmentationImage);
        rn -> set_point(i,j);
        nb = rn -> neighbourhood() & r;
        while (rn->point_count() < seed_size && nb.point_count()>0) {
		  *rn |= nb;
          nb = rn->neighbourhood() & r;
		} 
        if (rn->point_count() >= seed_size && dn < dsize && n < size) {
		  d[dn] = create(*rn, type);
          d[dn] -> normals = normals;
          if (d[dn]->error() < d[dn]->max_error()) {
            d[dn] -> print();
			segmentation::drawer->prepare(d[dn]->mmodel);
			segmentation::drawer->prepare(rn);
            handle[n++] = dn++;
		  } else { //failed
            delete d[dn];
		  }
		} else {}
        r &= (~(*rn));  
        delete rn;
	  }
	}
  qDebug() << "Placing seeds done";
  }

void segmentation::place_seeds(int seed_size, MODELTYPE type, segmentation *l)
{ region r(segmentationImage), *rn, nb(segmentationImage),rest(segmentationImage);
  int i,j,k;
  
  for (j = 0; j < segmentationImage->height(); j++)
    for (i = 0; i < segmentationImage->width(); i++)
      if (segmentationImage->valid_point(rest.get_point(i,j))) r.set_point(i,j);  

  rest.reset_all();

  for (k = 0; k < l -> n; k++)
    rest |= (*d[l->handle[k]]->mregion);  

  r &= ~rest;

    for (j = r.miny; j <= r.maxy; j++)
      for (i = r.minx; i <= r.maxx; i++) {
		if (!r.is_selected(i, j)) continue;
        if (r.included(i,j))
        { rn = new region(segmentationImage);
          rn -> set_point(i,j);
          nb = rn -> neighbourhood() & r;
          while (rn->point_count() < seed_size && nb.point_count()>0)
          { *rn |= nb;
            nb = rn->neighbourhood() & r;
            } 
          if (rn->point_count() >= seed_size && dn < dsize && n < size) 
          { d[dn] = create(*rn, type);
            d[dn] -> normals = normals;
            if (d[dn]->error() < d[dn]->max_error())
            { 
              d[dn] -> print();
			  segmentation::drawer->prepare(d[dn]->mmodel);
			  segmentation::drawer->prepare(rn);
              handle[n++] = dn++;
              } else
	    { //std::cout << "\nDescription creation failed - too big error " << d[dn]->error() << " ";
	      // std::cout << "\n";
              // d[dn]->mmodel->print();
              delete d[dn];
              }
            } else
	  { // std::cout << "\nDescription creation failed - not enough points for seed" << dn << ' ' << n
	    //   << ' ' << rn->point_count();
            }
          r &= (~(*rn));  
          delete rn;
          }
	  }
  }

void segmentation::import(MODELTYPE mtype, double *params) {
  region *rn = new region(segmentationImage);

  switch(mtype) {
  case CSQ:
	d[dn] = create(*rn, mtype);
	break;
  default:
	delete rn;
	return;
  }
  
  d[dn] -> normals = normals;
  d[dn]->mmodel->set_parameters(params);
  d[dn] -> print();
  
  //segmentation::drawer->prepare(d[dn]->mmodel);
  //segmentation::drawer->prepare(rn);
  handle[n++] = dn++;
  delete rn;
}

void segmentation::place_seeds(int seed_size, MODELTYPE type, segmentation *l, int maxseeds)
{ region r(segmentationImage), *rn, nb(segmentationImage),rest(segmentationImage);
  int i,j,k,hd = 0;
  
  r.reset_all();
  for (j = 0; j < segmentationImage->height(); j++)
    for (i = 0; i < segmentationImage->width(); i++)
      if (segmentationImage->valid_point(rest.get_point(i,j))) r.set_point(i,j);  

  rest.reset_all();

  for (k = 0; k < l -> n; k++)
    rest |= (*d[l->handle[k]]->mregion);  

  r &= ~rest;
  while (r.point_count() > 1.05 * seed_size * maxseeds)
  {  nb = r.boundary();
      r &= ~nb;
      }

    for (j = r.miny; j <= r.maxy; j++)
      for (i = r.minx; i <= r.maxx; i++)
        if (r.included(i,j) && hd < maxseeds)
        { rn = new region(segmentationImage);
          rn -> set_point(i,j);
          nb = rn -> neighbourhood() & r;
          while (rn->point_count() < seed_size && nb.point_count()>0)
          { *rn |= nb;
            nb = rn->neighbourhood() & r;
            } 
          if (rn->point_count() >= seed_size && dn < dsize && n < size) 
          { d[dn] = create(*rn, type);
            d[dn] -> normals = normals;
            if (d[dn]->error() < d[dn]->max_error())
            { 
              d[dn] -> print();
			  segmentation::drawer->prepare(d[dn]->mmodel);
			  segmentation::drawer->prepare(rn);
              handle[n++] = dn++;
              hd++;
              } else
	    { //std::cout << "\nDescription creation failed - too big error " << d[dn]->error() << " ";
	      // std::cout << "\n";
              // d[dn]->mmodel->print();
              delete d[dn];
              }
            } else
	  { // std::cout << "\nDescription creation failed - not enough points for seed" << dn << ' ' << n
	    //   << ' ' << rn->point_count();
            }
          r &= (~(*rn));  
          delete rn;
          }
  }

symatrix& segmentation::get_sel_matrix()
{ symatrix *m;
  int i,j;
  region r(segmentationImage);
  
  m = new symatrix(n);  
  for (i = 0; i < n; i++)
    for (j = i; j < n; j++)
      if (i == j)
      { m->el(i,i) = d[handle[i]] -> MSC();
        } else
      { m->el(i,j) = m->el(j,i) = d[handle[i]]->MSC(d[handle[j]]);
        }
  
  return(*m);
  }
  
void segmentation::selection()
{ 
  int i, mind, max_ind,ind,j;
  segmentation *in, *not_in;
  double max_delta_f, delta_f;
  symatrix sm(n);

  //std::cout << *this;
  if (n > 0)
  { throw_away();

    sm = get_sel_matrix();
    std::cout << "\n" << sm << "\n";
    in = new segmentation(size);
    not_in = new segmentation(size);
    for (i = 0; i < n; i++) 
    { not_in->handle[i] = i;            // not_in holds indices of i-th handlers not included
      in -> handle[i] = 0;              // only if in->handle[i]!=0 then handle[i] included
      }
    not_in->n = in->n = n;
    // std::cout << "\nInclude models:";
    do {
	  max_ind = not_in -> handle[mind = 0];
      max_delta_f = sm.el(max_ind,max_ind);
      for (i = 1; i < not_in->n; i++) {
		ind = not_in -> handle[i];
        delta_f = sm.el(ind,ind);
        if (delta_f > max_delta_f) {
		  max_delta_f = delta_f;
          max_ind = ind;
          mind = i;
		}
	  }
      if (max_delta_f > 0.0) {
		in->handle[max_ind] = 1;
        not_in->handle[mind] = not_in->handle[--not_in->n];
        for (j = 0; j < not_in->n; j++) {
		  ind = not_in->handle[j];
          sm.el(ind,ind) += 2*sm.el(ind,max_ind);
		}
        // std::cout << ' ' << max_ind << ' ' << max_delta_f;
        // std::cout.flush();  
	  }  
	} while (max_delta_f > 0.0 && not_in->n);
    for (i = 0; i < n; i++)
      if (!in->handle[i]) {
		delete d[handle[i]];
		d[handle[i]] = NULL;
	  }      
    throw_away();
    // std::cout << *this;
    // std::cout.flush();
    delete in;
    delete not_in;
    }    
  printf("\nList dn = %d",dn);
  fflush(stdout);
  }  

// void segmentation::draw(int mode)
// { int i;
//   region r(segmentationImage);
//   throw_away();

//   for (i = 0; i < n; i++)
//   { if (mode)
//     { glLoadName(handle[i]);
//       printf("Name %d loaded ",handle[i]);
//       }

//     d[handle[i]] -> draw();
//     r = d[handle[i]]->mregion->boundary();
//     r.draw(140);
//     }
  
//   glFlush();
//   }  
  
void segmentation::grow(int iterations)
{ int i,k;
  
  throw_away();
  for (k = 0; k < iterations; k++)
    for (i = 0; i < n; i++)
    { 
      if (d[handle[i]]->can_grow())
      { //printf("\nBEGINGROW");
	//pmem();
	
        d[handle[i]]->grow();
	  
        //pmem();
	//printf("\nENDGROW");
        
        //d[handle[i]]->draw();
        // pmem();
        // printf("\nENDDRAW");
	//        std::cout << "Description error: " << d[handle[i]]->error() << "\n";
	//        d[handle[i]]->mmodel->print();
        }
      }
  }  
    
description* segmentation::create(region& r, MODELTYPE type, model *m, int camera, double cX , double cY, double cZ)
{ double coord[3];
  description *d1;  

  d1 = 0;
  if (camera)
  { coord[0] = r.theImage->camOrientationX;
     coord[1] = r.theImage->camOrientationY;
     coord[2] = r.theImage->camOrientationZ;
     r.theImage->camOrientationX = cX;
     r.theImage->camOrientationY = cY;
     r.theImage->camOrientationZ = cZ;
     }

   switch (type)
  { case CPLANE: d1 = (new planar_d(r)); break;
    case CSQ: d1 = (new sq_d(r,m)); break;
    case CSURFACE2: d1 = (new surface2_d(r)); break;
    case CSPHERE: d1 = (new sphere_d(r)); break;
    case CCYLINDER: d1 = (new cylinder_d(r,m)); break;
    case CCONE: d1 = (new cone_d(r,m)); break;
    case CTORUS: d1 = (new torus_d(r,m)); break;
    case CASQ: d1 = (new asq_d(r,m)); break;
    case CTSQ: d1 = (new tsq_d(r,m)); break;
    case CBSQ: d1 = (new bsq_d(r,m)); break;
    }

   if (camera)
   {  r.theImage->camOrientationX = coord[0];
       r.theImage->camOrientationY = coord[1];
       r.theImage->camOrientationZ = coord[2];
       }

  return(d1);

  }  
    
void segmentation::set_sel_const(double kk2, double kk3)
{ k2 = kk2;
  k3 = kk3;
  }  
    
std::ostream& operator<<(std::ostream& s, segmentation& l)
{ int i;
  s << "\nList size " << l.n << "\nList elements:";
  for (i = 0; i < l.n; i++) s << ' ' << l.handle[i];
  return(s);
  }


double segmentation::calculateDeviation(int without)
{ int i,j,k,np = 0;
  double df,x = 0.0;

  for (k = 0; k < n; k++)
    if (k != without)
    { for (j = d[handle[k]]->mregion->miny; j <= d[handle[k]]->mregion->maxy; j++)
         for (i = d[handle[k]]->mregion->minx; i <= d[handle[k]]->mregion->maxx; i++)
            if (d[handle[k]]->mregion->included(i,j))
            { df = d[handle[k]]->mmodel->distance(d[handle[k]]->mregion->get_point(i,j));
              x += df * df;
              }
      np += d[handle[k]]->mregion->point_count();
      }
  return(x/(np-1));
  }

void segmentation::simpleSelection()
{ int a[n];
  int i,k = 0;
  double ds;

  gsSimple gs(this);
  gs.solve();
  ds = gs.solution(a);
  for (i = 0; i < n; i++)
     if (!a[i]) 
     { delete d[handle[i]];
       d[handle[i]] = NULL;
       } else k++;
  printf("\n%d out of %d simple descriptions selected",k,n); fflush(stdout);
  throw_away();
  }

void segmentation::simpleRAS(int seed_size, int fgi, int iter, MODELTYPE type)
{ int i,j,m;
  region r(segmentationImage);
  place_seeds(seed_size, type);
  m = fgi;
  for (i = 0; i < iter; i++)
  { for (j = 0; j < m; j++)
    { grow(1);
      printf("\nIteration %d/%d completed",j+1,m);
      fflush(stdout);
      }
    simpleSelection();
    m *= 2;
    }
  r.reset_all();
  for (i = 0; i < n; i++) r |= (*d[handle[i]]->mregion);
  printf("\n%d out of %d points included in simple RAS results",r.point_count(),
        segmentationImage->points());
  fflush(stdout);
  }
    
void segmentation::recover_and_select(int seed_size, int fgi, int iter, MODELTYPE type)
{ int i,n,k = 0;
  
  place_seeds(seed_size, type);
  n = fgi;

  for (i = 0; i < iter; i++)
	{ //segmentationImage->drawGL(); //TODO: migrate
	  //glFlush();
	  grow(n);
	  n *= 2;
    // std::cout << "\nSelecting...";
    // std::cout.flush();
    // segmentationImage->drawGL(); // TODO: migrate
    selection();
    // if (showtype != RSGEOMVIEW) print_usage();
    }

  }  
  
segmentation::segmentation(FILE *f, image *im_1, image *norm_1)
{ int i;
  MODELTYPE model_type;
  int dwrong, camera_set = 0;

  segmentationImage = im_1;
  normals = norm_1;

  if (d == NULL)
  { fread(&dsize, sizeof(int), 1, f);
     if (dsize < 0)
     {  fread(&dsize, sizeof(int), 1, f);
         camera_set = 1;
         } camera_set = 0;

    d = new description*[dsize];
    dn = 0;
    } else fread(&dwrong,sizeof(int),1,f);

  fread(&size, sizeof(int), 1, f);
  fread(&n, sizeof(int), 1, f);
  fread(&k2, sizeof(double), 1, f);
  fread(&k3, sizeof(double), 1, f);
  fread(&(model::minx),sizeof(double),1,f);
  fread(&(model::miny),sizeof(double),1,f);
  fread(&(model::dx),sizeof(double),1,f);
  fread(&(model::dy),sizeof(double),1,f);
  fread(&(model::image_width),sizeof(int),1,f);
  fread(&(model::image_height),sizeof(int),1,f);
  handle = new int[size];
  for (i = 0; i < n; i++) 
  { handle[i] = dn + i;
    fread(&model_type,sizeof(MODELTYPE),1,f);
    switch(model_type)
    { case CPLANE: d[handle[i]] = new planar_d(f,segmentationImage,normals, camera_set);
                   break;
      case CSQ: d[handle[i]] = new sq_d(f,segmentationImage,normals, camera_set);
                break;
      case CSURFACE2: d[handle[i]] = new surface2_d(f,segmentationImage,normals, camera_set);
                      break;
      case CSPHERE: d[handle[i]] = new sphere_d(f,segmentationImage,normals, camera_set);
                    break;
      case CCYLINDER: d[handle[i]] = new cylinder_d(f,segmentationImage,normals, camera_set);
                      break;
      case CCONE: d[handle[i]] = new cone_d(f,segmentationImage,normals, camera_set);
                      break;
      case CTORUS: d[handle[i]] = new torus_d(f,segmentationImage,normals, camera_set);
                      break;
	  case CASQ: d[handle[i]] = new asq_d(f,segmentationImage,normals, camera_set);
                      break;
	  case CTSQ: d[handle[i]] = new tsq_d(f,segmentationImage,normals, camera_set);
		              break;
	  case CBSQ: d[handle[i]] = new bsq_d(f,segmentationImage,normals, camera_set);
                      break;
      }
    d[handle[i]] -> normals = normals;
    }
  dn += n;  
  usingList++;
  }  

void segmentation::load(FILE *f, image *im_1, image *norm_1)
{ int i;
  MODELTYPE model_type;
  int dwrong, camera_set = 0;

  segmentationImage = im_1;
  normals = norm_1;

  if (d == NULL)
  { fread(&dsize, sizeof(int), 1, f);
     if (dsize < 0)
    {  fread(&dsize,sizeof(int),1,f);
        camera_set = 1;
        } else camera_set = 0;

    d = new description*[dsize];
    dn = 0;
    } else 
     {  fread(&dwrong,sizeof(int),1,f);
         if (dwrong < 0) 
         {  fread(&dwrong, sizeof(int), 1, f);
             camera_set = 1;
             } else camera_set = 0;
        }

  fread(&size, sizeof(int), 1, f);
  fread(&n, sizeof(int), 1, f);
  fread(&k2, sizeof(double), 1, f);
  fread(&k3, sizeof(double), 1, f);
  fread(&(model::minx),sizeof(double),1,f);
  fread(&(model::miny),sizeof(double),1,f);
  fread(&(model::dx),sizeof(double),1,f);
  fread(&(model::dy),sizeof(double),1,f);
  fread(&(model::image_width),sizeof(int),1,f);
  fread(&(model::image_height),sizeof(int),1,f);
  delete [] handle;
  handle = new int[size];
  for (i = 0; i < n; i++) 
  { handle[i] = dn + i;
    fread(&model_type,sizeof(MODELTYPE),1,f);
    switch(model_type)
    { case CPLANE: d[handle[i]] = new planar_d(f,segmentationImage,normals, camera_set);
                   break;
      case CSQ: d[handle[i]] = new sq_d(f,segmentationImage,normals, camera_set);
                break;
      case CSURFACE2: d[handle[i]] = new surface2_d(f,segmentationImage,normals, camera_set);
                      break;
      case CSPHERE: d[handle[i]] = new sphere_d(f,segmentationImage,normals, camera_set);
                    break;
      case CCYLINDER: d[handle[i]] = new cylinder_d(f,segmentationImage,normals, camera_set);
                      break;
      case CCONE: d[handle[i]] = new cone_d(f,segmentationImage,normals, camera_set);
                      break;
      case CTORUS: d[handle[i]] = new torus_d(f,segmentationImage,normals, camera_set);
                      break;
      case CASQ: d[handle[i]] = new asq_d(f,segmentationImage,normals, camera_set);
                      break;
	  case CTSQ: d[handle[i]] = new tsq_d(f,segmentationImage,normals, camera_set);
                      break;
	  case CBSQ: d[handle[i]] = new bsq_d(f,segmentationImage,normals, camera_set);
                      break;
      }
    d[handle[i]] -> normals = normals;
    }
  dn += n;  
  usingList++;
  }  

// Make sure you call the next function after constructor
// and before doing anything with segmentation

void segmentation::init_list(image *im, image *norm)
{ segmentationImage = im;
  normals = norm;
}                                                     

void segmentation::delete_description(int i)
{ delete d[handle[i]];
  d[handle[i]] = NULL;
  throw_away();
  }

void segmentation::adjacency_graph(int **graph, int limit)
{ int i,j,k;
  region r(segmentationImage),br(segmentationImage),r1(segmentationImage);


  throw_away();
  *graph = new int[n*n];  

  for (i = 0; i < n; i++)
  { r = *d[handle[i]]->mregion;
    r |= r.neighbourhood();
    (*graph)[i*(n+1)] = 0;  
    for (j = i + 1; j < n; j++)
    { r1 = *d[handle[j]]->mregion;
      r1 |= r1.neighbourhood();
      br = r & r1;
      if (br.point_count() > limit) k = 1;
        else k = 0;
      (*graph)[i+j*n] = (*graph)[j+i*n] = k;
      }  
    }
  }

void segmentation::exclude_intersections(double ep)
{ int i,j,c1,c2,cint;
  region r(segmentationImage),rn(segmentationImage),isc(segmentationImage);
  MODELTYPE m;
  int *to_delete;
  description *dtemp;

  to_delete = new int[n];
  for (i = 0; i < n; i++) to_delete[i] = 0;

  isc.reset_all();
  for (i = 0; i < n; i++)
    for (j = i + 1; j < n; j++)
    { r = *d[handle[i]]->mregion & *d[handle[j]]->mregion;
      cint = r.point_count();
      c1 = d[handle[i]]->mregion->point_count();
      c2 = d[handle[j]]->mregion->point_count();
      if (cint > ep*c1 && c1 < c2)
      { to_delete[i] = 1;
        printf("\nSimilar description %d deleted: %d %d %d",i,c1,c2,cint);
	} else if (cint > ep*c2 && c2 <= c1)
	  { to_delete[j] = 1; 
            printf("\nSimilar description %d deleted: %d %d %d",j,c1,c2,cint);
	    } else
              { isc |= r;
	        }
      }

  rn = ~isc;

  for (i = 0; i < n; i++)
    if (to_delete[i])
    { delete d[handle[i]];
      d[handle[i]] = NULL;
      } else
      { *d[handle[i]]->mregion &= rn;
	}

  throw_away();

  for (i = 0; i < n; i++) 
  { r = *d[handle[i]]->mregion;
    m = d[handle[i]]->mmodel->what_model();
    dtemp = create(r,m,d[handle[i]]->mmodel);
    delete d[handle[i]];
    d[handle[i]] = dtemp;
    d[handle[i]] -> normals = normals;
    }
  //draw();
  }

int segmentation::join_compatible(int constr)
{ int i,j,k,n2,ind,stat = 0,i1,i2,tok,stat1 = 0;
  struct point p;
  region **r;
  double d1,d2,d3,derr,dmin,dist;

  r = new region*[n];

  for (i = 0; i < n; i++) 
  { r[i] = new region(segmentationImage); 
    r[i]->reset_all();
    }
  for (j = 0; j < segmentationImage->points(); j++)
  { p = segmentationImage -> pixel(j);
    n2 = (p.n >> 1);
    for (k = 0; k < n2; k++)
    { i1 = p.neigh[2*k];
      i2 = p.neigh[2*k+1];
      if (j < i1 && j < i2 && triangle_free(j,i1,i2))
      { if (constr) 
        { dmin = constr * pp_max_err;
          dist = constr * pp_max_dist;
          } else 
        { dmin = 1e10;  // infinity
          dist = 1e10; // infinity
          }
        ind = -1;
        tok = 0;
        for (i = 0; i < n; i++)
        { if (d[handle[i]]->mregion->triangle_near(j,i1,i2))
          { d1 = d[handle[i]]->mmodel->distance(p);
	    d2 = d[handle[i]]->mmodel->distance(segmentationImage->pixel(i1));
            d3 = d[handle[i]]->mmodel->distance(segmentationImage->pixel(i2));
            derr = (d1 + d2 + d3) / 3.0;
            if (d1 < dist && d2 < dist && d3 < dist && derr < dmin)
	    { dmin = derr;
              ind = i;
	      }
            tok = 1;
	    }
	  }        
        if (ind > -1)
        { r[ind]->set_point(j);
          r[ind]->set_point(i1);
          r[ind]->set_point(i2);
          stat++;
	  }
        if (tok) stat1++;
        }
      }      
    }
  
  for (i = 0; i < n; i++)
  { *d[handle[i]]->mregion |= *r[i];
    delete r[i];
    }
  delete [] r;
  printf("\n%d triangles out of %d scheduled in postprocessing",stat,stat1); 
  return(stat);
  }
  
int segmentation::triangle_free(int i1, int i2, int i3)
{ int i;
  region *r;
  for (i = 0; i < n; i++)
  { r = d[handle[i]] -> mregion;
    if (r->included(i1) && r->included(i2) && r->included(i3)) return(0);
    }
  return(1);
  }

void segmentation::label_points(FILE *f)
{ int i,j,np,k;
  int *in;

  in = new int[n];
  np = segmentationImage -> points();
  
  for (i = 0; i < np; i++)
  { k = 0;
    for (j = 0; j < n; j++)
      if (d[handle[j]]->mregion->included(i)) in[k++] = j+1;
    fprintf(f,"%d ",k);
    for (j = 0; j < k; j++) fprintf(f,"%d ",in[j]);
    if (i < np - 1) fprintf(f,"\n");        
    }

  delete [] in;
  }

void segmentation::fprint(FILE *f)
{ int i;

  for (i = 0; i < n; i++)
  { d[handle[i]]->mmodel->fprint(f);
    if (i < n - 1) fprintf(f,"\n");
    }
  }

void segmentation::model_register(int *k, segmentation &l, symatrix &A)
{ symatrix M(3),M1(3);
  vector a(3),b(3),c(3),res(3),dd(3);
  int i;

  A.el(3,3) = 1;
  for (i = 0; i < 3; i++) A.el(3,i) = 0;

  for (i = 0; i < 3; i++) 
    if (d[handle[k[i]]]->mmodel->what_model() != CPLANE || 
      l.d[l.handle[k[i+3]]]->mmodel->what_model() != CPLANE)
        fprintf(stderr,"\nNot correct type of model for registration!");
   
  for (i = 0; i < 3; i++)
  { d[handle[k[i]]]->mmodel->get_normal_vector(M.el(i,0),M.el(i,1),M.el(i,2));
    l.d[l.handle[k[i+3]]]->mmodel->get_normal_vector(a.el(i),b.el(i),c.el(i));
    }

  M1 = M.inverse(i);

  res = M1 * a;
  for (i = 0; i < 3; i++) A.el(i,0) = res.el(i);
  res = M1 * b;
  for (i = 0; i < 3; i++) A.el(i,1) = res.el(i);
  res = M1 * c;
  for (i = 0; i < 3; i++) A.el(i,2) = res.el(i);

  std::cout << "\n"  << a << "\n" << b << "\n" << c;
  
  for (i = 0; i < 3; i++) 
    dd.el(i) = l.d[l.handle[k[i+3]]]->mmodel->model_distance()-d[handle[k[i]]]->mmodel->model_distance();

  res = M1 * dd;
  for (i = 0; i < 3; i++) A.el(i,3) = res.el(i);

  // std::cout << "\n" << M.determinant() << " " << M1.determinant() << "\n"; 
  // std::cout << M << "\n" << M1 << "\n" << M * M1 << "\n" << A;

  std::cout << M1;
  }

void segmentation::image_transform(symatrix &A)
{ struct point p;
  vector v(4),r(4),vn(3),rn(3);
  vector3D tr;
  symatrix C(4),B(3);
  FILE *f;
  double norm;

  int i,j,np = segmentationImage->points();

  C.identity();

  for (i = 0; i < 3; i++)
  { for (j = 0; j < 3; j++) 
    { B.el(i,j) = A.el(i,j);
      C.el(i,j) = A.el(i,j);
      }
    tr.el(i) = A.el(i,3);
    }

  v.el(3) = r.el(3) = 1.0;

  f = fopen("xxx.vtx","w");

  for (i = 0; i < np; i++)
  { p = segmentationImage->pixel(i);
    v.el(0) = p.x;
    v.el(1) = p.y;
    v.el(2) = p.z;
    r = A * v;
    segmentationImage->set_pixel(i,0,r.el(0),r.el(1),r.el(2));
    fprintf(f,"%lf %lf %lf\n",r.el(0),r.el(1),r.el(2));
    if (normals != NULL)
    { vn.el(0) = p.nx;
      vn.el(1) = p.ny;
      vn.el(2) = p.nz;
      rn = B * vn;
      norm = rn.norm();
	if (norm>0.9 && norm < 1.0) 
     {  rn.el(0) /= norm;
         rn.el(1) /= norm;
         rn.el(2) /= norm;
         } 
      segmentationImage->set_normal(i,rn.el(0),rn.el(1),rn.el(2));
      normals->set_pixel(i,0,rn.el(0),rn.el(1),rn.el(2));
      }
    } 
  
  fclose(f);
  segmentationImage->extremes();

  for (i = 0; i < n; i++)
  { d[handle[i]]->mmodel->rotate(C);
    d[handle[i]]->mmodel->translate(tr);
    }

  std::cout << "Transformation: " << B << "\n" << tr;
  }

void segmentation::rotate(symatrix &A)
{ struct point p;
  vector v(4),r(4),vn(3),rn(3), cam1(4), cam2(4);
  symatrix B(3);
  double norm;

  int i,j,np = segmentationImage->points();

  cam1.el(0) = segmentationImage->camOrientationX;
  cam1.el(1) = segmentationImage->camOrientationY;
  cam1.el(2) = segmentationImage->camOrientationZ;
  cam1.el(3) = 1.0;
  cam2 = A * cam1;
  segmentationImage->camOrientationX = cam2.el(0);
  segmentationImage->camOrientationY = cam2.el(1);
  segmentationImage->camOrientationZ = cam2.el(2);

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++) B.el(i,j) = A.el(i,j);

  v.el(3) = r.el(3) = 1.0;

  for (i = 0; i < np; i++)
  { p = segmentationImage->pixel(i);
    v.el(0) = p.x;
    v.el(1) = p.y;
    v.el(2) = p.z;
    r = A * v;
    segmentationImage->set_pixel(i,0,r.el(0),r.el(1),r.el(2));
    if (normals != NULL)
    { vn.el(0) = p.nx;
      vn.el(1) = p.ny;
      vn.el(2) = p.nz;
      rn = B * vn;
      norm = rn.norm();
	if (norm>0.9 && norm < 1.0) 
     {  rn.el(0) /= norm;
         rn.el(1) /= norm;
         rn.el(2) /= norm;
         } 
      segmentationImage->set_normal(i,rn.el(0),rn.el(1),rn.el(2));
      normals->set_pixel(i,0,rn.el(0),rn.el(1),rn.el(2));
      }
    } 
  
  for (i = 0; i < n; i++)
  { d[handle[i]]->mmodel->rotate(A);
    }

  segmentationImage->extremes();
  }

void segmentation::translate(vector3D& v)
{ int i,np = segmentationImage->points();
  struct point p;

  for (i = 0; i < np; i++)
  { p = segmentationImage->pixel(i);
    p.x += v.el(0);
    p.y += v.el(1);
    p.z += v.el(2);
    segmentationImage->set_pixel(i,0,p.x,p.y,p.z);
    }

  for (i = 0; i < n; i++)
  { d[handle[i]]->mmodel->translate(v);
    }

  segmentationImage->extremes();

  }

void segmentation::ICP(double dist, segmentation &l, symatrix &A,region &r1, region &r2)
{ struct point po,pt;
  int i,j,mi,n1,n2,nt = 0;
  double *vo,*vt,dm,d;
  symatrix B(3);
  vector3D vr;

  n1 = segmentationImage->points();
  n2 = l.segmentationImage->points();
  vo = new double[4*n1];
  vt = new double[4*n1];
  r1.reset_all();
  r2.reset_all();

  for (i = 0; i < n1; i++)
  { po = segmentationImage->pixel(i);
    pt = l.segmentationImage->pixel(0);
    mi = 0;
    dm = distance(po,pt);
    for (j = 1; j < n2; j++)
    { pt = l.segmentationImage->pixel(j);
      d = distance(po,pt);
      if (d < dm)
      { dm = d;
        mi = j;
	}
      }
    if (dm < dist)
    { vo[3*nt] = po.x; 
      vo[3*nt+1] = po.y;
      vo[3*nt+2] = po.z;
      pt = l.segmentationImage->pixel(mi);
      vt[3*nt] = pt.x;
      vt[3*nt+1] = pt.y;
      vt[3*nt+2] = pt.z; 
      nt++;
      r1.set_point(i);
      r2.set_point(mi);
      }   
    }
  
  printf("\n%d %d comparable points\n",nt,r2.point_count());
  p_rotation(vt,vo,nt,A);
  for (i = 0; i < 3; i++)
     for (j = 0; j < 3; j++)
        B.el(i,j) = A.el(i,j);

  for (i = 0; i < nt; i++)
  {  vr.el(0) = vt[3*i];
      vr.el(1) = vt[3*i+1];
      vr.el(2) = vt[3*i+2];
      vr = B * vr;
      vt[3*i] = vr.el(0);
      vt[3*i+1] = vr.el(1);
      vt[3*i+2] = vr.el(2);
      }
  translation(vt,vo,nt,vr);

  A.set_submatrix(vr,0,3);  

  delete [] vo;
  delete [] vt;  
  }

void segmentation::ICRP(double dist, segmentation &l, symatrix &A,region &r1, region &r2)
{ struct point po,pt,pr;
  int i,j,k,mi,n1,n2,nt = 0, mr;
  double *vo,*vt,dm,d, dr;
  symatrix B(3);
  vector3D vr;

  n1 = segmentationImage->points();
  n2 = l.segmentationImage->points();
  vo = new double[4*n1];
  vt = new double[4*n1];
  r1.reset_all();
  r2.reset_all();

  for (i = 0; i < n1; i++)
  { po = segmentationImage->pixel(i);
    pt = l.segmentationImage->pixel(0);
    mi = 0;
    dm = distance(po,pt);
    for (j = 1; j < n2; j++)
    { pt = l.segmentationImage->pixel(j);
       d = distance(po,pt);
       if (d < dm)
       { dm = d;
          mi = j;
	    }
       }

     pt = l.segmentationImage->pixel(mi);
     // reciprocal check
    pr = segmentationImage->pixel(0);
    mr = 0;
    dr = distance(pr,pt);
    for (k = 1; k < n1; k++)
    {  pr = segmentationImage->pixel(k);
        d = distance(pr,pt);
        if (d < dr)
        {  dr = d;
            mr = k;
            }
        }

    if (dm < dist && mr == i)
    { vo[3*nt] = po.x; 
      vo[3*nt+1] = po.y;
      vo[3*nt+2] = po.z;
      pt = l.segmentationImage->pixel(mi);
      vt[3*nt] = pt.x;
      vt[3*nt+1] = pt.y;
      vt[3*nt+2] = pt.z; 
      nt++;
      r1.set_point(i);
      r2.set_point(mi);
      }   
    }
  
  printf("\n%d %d comparable points\n",nt,r2.point_count());
  p_rotation(vt,vo,nt,A);
  for (i = 0; i < 3; i++)
     for (j = 0; j < 3; j++)
        B.el(i,j) = A.el(i,j);

  for (i = 0; i < nt; i++)
  {  vr.el(0) = vt[3*i];
      vr.el(1) = vt[3*i+1];
      vr.el(2) = vt[3*i+2];
      vr = B * vr;
      vt[3*i] = vr.el(0);
      vt[3*i+1] = vr.el(1);
      vt[3*i+2] = vr.el(2);
      }
  translation(vt,vo,nt,vr);

  A.set_submatrix(vr,0,3);  

  delete [] vo;
  delete [] vt;  
  }

void segmentation::ProjectICP(double dist, int nk, int *k1, int *k2, segmentation &l,symatrix& A,region& r2)
{ double *vo,*vt,mindist,dist1, maxdist;
  int n1,n2,ind,i,j,nt;
  struct point p,p1;
  vector3D vr;
  symatrix B(3);

  n1 = segmentationImage->points();
  n2 = l.segmentationImage->points();
  vo = new double[3*(n1+n2)];
  vt = new double[3*(n1+n2)];
  nt = 0;  
  maxdist = 0.0;

  for (i = 0; i < n1; i++)
  { mindist = dist;
     ind = -1;
     p = segmentationImage->pixel(i);
     for (j = 0; j < nk; j++)
        if (d[handle[k1[j]]]->mregion->included(i))            
        { dist1 = l.d[l.handle[k2[j]]]->mmodel->distance(p);  
           if (dist1 < mindist)
           {  mindist = dist1;
              ind = j;
              }
           }   

    if (ind > -1)     
    { p1 = l.d[l.handle[k2[ind]]]->mmodel->project(p);    
      dist1 = distance(p,p1); 
      if (dist1 > maxdist) maxdist = dist1;
      if (dist1 < dist)
      {  vo[3*nt] = p.x;
         vo[3*nt+1] = p.y;
         vo[3*nt+2] = p.z;
         vt[3*nt] = p1.x;
         vt[3*nt+1] = p1.y;
         vt[3*nt+2] = p1.z;
         nt++;
         r2.set_point(i);
         }
      }
    }

//  reverse

  for (i = 0; i < n2; i++)
  { mindist = dist;
     ind = -1;
     p = l.segmentationImage->pixel(i);
     for (j = 0; j < nk; j++)
        if (l.d[l.handle[k2[j]]]->mregion->included(i))            
        { dist1 = d[handle[k1[j]]]->mmodel->distance(p);  
           if (dist1 < mindist)
           {  mindist = dist1;
              ind = j;
              }
           }   

    if (ind > -1)     
    { p1 = d[handle[k1[ind]]]->mmodel->project(p);    
      dist1 = distance(p,p1); 
      if (dist1 > maxdist) maxdist = dist1;
      if (dist1 < dist)
      {  vo[3*nt] = p1.x;
         vo[3*nt+1] = p1.y;
         vo[3*nt+2] = p1.z;
         vt[3*nt] = p.x;
         vt[3*nt+1] = p.y;
         vt[3*nt+2] = p.z;
         nt++;
         //r2.set_point(i);
         }
      }
    }

  printf("\n %d comparable points\n",nt);
  printf("Max distance:%lf\n\n", maxdist);
  p_rotation(vt,vo,nt,A);
  for (i = 0; i < 3; i++)
     for (j = 0; j < 3; j++)
        B.el(i,j) = A.el(i,j);

  for (i = 0; i < nt; i++)
  {  vr.el(0) = vt[3*i];
      vr.el(1) = vt[3*i+1];
      vr.el(2) = vt[3*i+2];
      vr = B * vr;
      vt[3*i] = vr.el(0);
      vt[3*i+1] = vr.el(1);
      vt[3*i+2] = vr.el(2);
      }
  
   translation(vt,vo,nt,vr);
  A.set_submatrix(vr,0,3);  

  delete [] vo;
  delete [] vt;

  }

void segmentation::MedioniICP(double dist,segmentation &l,symatrix& A, region &r1, region& r2)
{ struct point po,pt;
  int i,j,mi,n1,n2,nt = 0;
  double *vo,*vt,dm,d;
  vector3D vr;
  symatrix B(3);
  plane *pl;

  n1 = segmentationImage->points();
  n2 = l.segmentationImage->points();
  vo = new double[4*n1];
  vt = new double[4*n1];
  r1.reset_all();
  r2.reset_all();

  for (i = 0; i < n1; i++)
  { po = segmentationImage->pixel(i);
    pt = l.segmentationImage->pixel(0);
    mi = 0;
    dm = distance(po,pt);
    for (j = 1; j < n2; j++)
    { pt = l.segmentationImage->pixel(j);
      d = distance(po,pt);
      if (d < dm)
      { dm = d;
         mi = j;
  	  }
     }
    if (dm < dist)
    { vo[3*nt] = po.x; 
      vo[3*nt+1] = po.y;
      vo[3*nt+2] = po.z;
      pt = l.segmentationImage->pixel(mi);
	d = sqrt(pt.nx * pt.nx + pt.ny * pt.ny + pt.nz * pt.nz);
	if (d < 0.5) 
      { pt.nx = -l.segmentationImage->camOrientationX;
         pt.ny = -l.segmentationImage->camOrientationY;
         pt.nz = -l.segmentationImage->camOrientationZ;

	   }
      pl = new plane(pt.nx,pt.ny,pt.nz, -pt.x * pt.nx - pt.ny * pt.y - pt.nz * pt.z);
	pt = pl->project(po);
	
	//printf(" %lf",distance(pt,po));
	delete pl;
      vt[3*nt] = pt.x;
      vt[3*nt+1] = pt.y;
      vt[3*nt+2] = pt.z; 
      nt++;
      r1.set_point(i);
      r2.set_point(mi);
      }   
    }
  
  printf("\n%d %d comparable points\n",nt,r2.point_count());
  
  p_rotation(vt,vo,nt,A);
  for (i = 0; i < 3; i++)
     for (j = 0; j < 3; j++)
        B.el(i,j) = A.el(i,j);

  for (i = 0; i < nt; i++)
  {  vr.el(0) = vt[3*i];
      vr.el(1) = vt[3*i+1];
      vr.el(2) = vt[3*i+2];
      vr = B * vr;
      vt[3*i] = vr.el(0);
      vt[3*i+1] = vr.el(1);
      vt[3*i+2] = vr.el(2);
      }
  
   translation(vt,vo,nt,vr);
  A.set_submatrix(vr,0,3);  

  delete [] vo;
  delete [] vt;  
  }

void segmentation::delete_wrong()
{ int i;

  for (i = 0; i < n; i++)
    if (!d[handle[i]]->mmodel->model_ok())
    { delete d[handle[i]];
      d[handle[i]] = NULL;
      }
  throw_away();
  }

double segmentation::goodness()
{ int i,j;
  double sum = 0.0,q;
  
  for (i = 0; i < n; i++)
  { //if (n < 15) printf("\n\n");
     for (j = 0; j < n; j++)
        if (d[handle[i]] != NULL && d[handle[j]] != NULL) 
        { q = d[handle[i]]->MSC(d[handle[j]]);
           sum += q;
           //if (n < 15) printf("%lf ",q);
           }
     }
  if (n < 15) 
  { //printf("\n");
     //fflush(stdout);
     }
  return(sum);
  }

double segmentation::goodnessOfOne(int k)
{ int i;
  double sum = 0.0;
  
  for (i = 0; i < n; i++)
     if (i == k) sum += d[handle[k]]->MSC();
       else sum += 2*d[handle[k]]->MSC(d[handle[i]]);
  return(sum);
  }

double segmentation::goodnessOfTwo(description *d1, description *d2, int k, int l)
{ int i;
  double sum;
  description *d3;

  if (d1 == NULL)
  { d3 = d1;
    d1 = d2;
    d2 = d3;
    }
  if (d2 != NULL) sum = d1->MSC() + d2->MSC() + 2 * d1->MSC(d2);
    else sum = d1->MSC();

  for (i = 0; i < n; i++)
    if (i != k && i != l)
    { sum += 2 * d1->MSC(d[handle[i]]);
      if (d2 != NULL) sum += 2 * d2->MSC(d[handle[i]]);
      }
  return(sum);
  }

double segmentation::goodnessOfOne(description *d1, int k, int l)
{ return(goodnessOfTwo(d1,NULL,k,l));
  }

inline double segmentation::intersectionPenalty(int i, int j)
{ return(2*d[handle[i]]->MSC(d[handle[j]]));
  }

void segmentation::reset()
{ throw_away();
  n = 0;
  }

segmentation::~segmentation()
{ int i;
  usingList--;
  if (!usingList) 
  { for (i = 0; i < dn; i++) 
      if (d[i] != NULL) delete d[i];
        
    delete [] d;
    d = NULL;
    }
  delete [] handle;
  }
  

gsSimple::gsSimple(segmentation *l1) : gsBQP(l1->n)
{ l = l1;
  }

double gsSimple::el_source(int i, int j)
{ region r = *l->d[l->handle[i]]->mregion;

  if (i != j) r &= *l->d[l->handle[j]]->mregion;
  if (i == j) return(1.0*r.point_count());
    else return(-0.5*r.point_count());
  }

description** segmentation::d = NULL;
int segmentation::usingList = 0;  
image *normal = NULL;
double segmentation::pp_max_err = 0.0;
double segmentation::pp_max_dist = 0.0;
int segmentation::dsize = 0;
int segmentation::dn = 0;

