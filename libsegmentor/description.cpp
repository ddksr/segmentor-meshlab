// description.C

#include "description.h"
#include "chiMSC3D.h"

int useStatistics = 0;

overflowManager oM;


overflowManager::overflowManager(const char *name, int nmax)
{ int i,ix,iy;
   FILE *f;

   sprintf(filename,"%s",name);
   no = nmax;
   over = new int[no];
   for (i = 0; i < no; i++) over[i] = -1;

   f = fopen(filename,"rt");
   if (f != NULL)
   { while (fscanf(f,"%d %d",&ix,&iy) == 2) over[ix] = iy;
      fclose(f);
      }
   
   }

overflowManager::~overflowManager()
{ delete [] over;
   }

double overflowManager::P(int m, int n)
{ double rez = 1.0;
   int i,j,k;

   i = j = k = 1;
   while (i <= n || j <=m || k <= n-m)
   { if (i <= n) 
      { rez *= (i/2.0);
         i++;
         }

      if (rez >= 1e-5 || i > n)
     { if (j <= m)
        { rez /= j;
           j++;
           }
        if (k <= n-m)
        { rez /= k;
           k++;
           }
        }      
      }
 
   if (2*m == n) return(rez);
     else return(2*rez);
   }

int overflowManager::overflowP(int n, double thr)
{ int m = n / 2 + n % 2;
   double rez = 0.0;
   FILE *f;

   if (n < no && n > 0)
      if (over[n] != -1) return(over[n]);

   if (n < 2) return(n+1);

   while (rez <= thr)
   { rez += P(m,n);
      m++;
      }

  m--;
  if (n < no) over[n] = m;
  f = fopen(filename,"at");
  fprintf(f,"%d %d\n",n,m); 
  fclose(f);
  return(m);
   }


description::description(region& r)
{ camera_set = 1;
   cameraCoordinates[0] =r.theImage->camOrientationX;
   cameraCoordinates[1] =r.theImage->camOrientationY;
   cameraCoordinates[2] =r.theImage->camOrientationZ;
   can_grow_var = 1;
   allowed = new region(r.theImage);
   allowed->set_all();
   }

void description::grow(region *ra)
{ region r = mregion -> neighbourhood() & (*allowed);
  region rn = r;
  int i,j;
  model *from_old, *plain;
  double d1,d2;
  FILE *f;  
  char str[1000];
  double coord[3];

  if (can_grow() || ra != NULL)
  { if (ra == NULL)
    { for (j = r.miny; j <= r.maxy; j++)
        for (i = r.minx; i <= r.maxx; i++)
          if (r.included(i,j) && mregion->theImage->valid_point(mregion->get_point(i,j)))
            if (!compatible(i,j))
              r.reset_point(i,j); 
                else
              rn.reset_point(i,j);
       } else
       { r = *ra;
         rn.reset_all();
         }
    //r.draw(141);
    //rn.draw(140);        
    if (r.point_count())
    {  if (camera_set)
       { coord[0] = r.theImage->camOrientationX;
          coord[1] = r.theImage->camOrientationY;
          coord[2] = r.theImage->camOrientationZ;
          r.theImage->camOrientationX = cameraCoordinates[0];
          r.theImage->camOrientationY = cameraCoordinates[1];
          r.theImage->camOrientationZ = cameraCoordinates[2];
          }

      plain = mmodel -> improve(*mregion, r);
      from_old = mmodel -> improve(mmodel, *mregion, r);

      if (camera_set)
     {  r.theImage->camOrientationX = coord[0];
         r.theImage->camOrientationY = coord[1];
         r.theImage->camOrientationZ = coord[2];
         }

      d1 = error(plain);
      d2 = error(from_old);
      if (d1 < max_error() || d2 < max_error())
      { delete mmodel;
        if (d1 < d2) 
        { mmodel = plain;
	// cout << "\nFrom plain\n";
          delete from_old;
          } else 
        { mmodel = from_old;
	// cout << "\nFrom old\n";
          delete plain;
          }
        *mregion |= r;  
         if (useStatistics) analyze();                // check for statistics support
        } else
      { can_grow_var = 0;  // too big model error
        std::cout << "\nToo big error " << d1 << " " << d2 << " " << max_error() << "\n";
        f = fopen("segmentation.log","at");
        fprintf(f,"Too big error\n");
        mmodel->describe_model(str);
        fprintf(f,"Previous model: %s error=%lf\n",str,error(mmodel));
        plain->describe_model(str);
        fprintf(f,"From plain: %s error=%lf\n",str,d1);
        from_old->describe_model(str);
        fprintf(f,"From old: %s error=%lf\n",str,d2);       
        fclose(f);
        delete plain;
        delete from_old;
        }
      } else
    { can_grow_var = 0;   // no more points
      std::cout << "\nNo more compatible points";
      f = fopen("segmentation.log","at");
      fprintf(f,"No more compatible points\n");
      mmodel->describe_model(str);
      fprintf(f,"Model: %s error=%lf\n",str,error(mmodel));
      fclose(f);      
      }
    }   
  }

double description::total_error(region& r, int& card)
{ int i,j;
  double sum = 0.0;
  
  card = 0;
  for (j = r.miny; j <= r.maxy; j++)
    for (i = r.minx; i <= r.maxx; i++)
      if (r.included(i,j))
      { sum += mmodel->distance(r.get_point(i,j));
        card++;
        }
  return(sum);
        
  }  
  
double description::error(model *m)
{ int i,j,count = 0;
  double sum = 0;
  
  for (j = mregion -> miny; j <= mregion -> maxy; j++)
    for (i = mregion -> minx; i <= mregion -> maxx; i++)
      if (mregion -> included(i,j))
      { sum += m -> distance(mregion -> get_point(i,j));
        count++;
        }
  
  if (count) return(sum/count);
    else return(0.0);      
  }  
  
double description::noise_variance()
{ int i,j;
   double aver1 = 0.0, aver2 = 0.0, dis, rez;
   int card = 0;   
   
  for (j = mregion -> miny; j <= mregion -> maxy; j++)
    for (i = mregion -> minx; i <= mregion -> maxx; i++)
      if (mregion -> included(i,j))
      { dis = mmodel->signed_distance(mregion->get_point(i,j));
         aver1 += dis;
         aver2 += dis*dis;
         card++;
         }   

   printf("\nAverage: %lf",aver1/card);
   fflush(stdout);
   rez = aver2 / card - (aver1 / card) * (aver1 / card);
   return(rez);

   }

// void description::draw()
// { int i,j;
 
//   if (drawRegions) mregion -> draw();
//   mmodel -> draw();
//   if (drawRegions)
//   { glBegin(GL_POINTS);
//      for (j = mregion -> miny; j <= mregion -> maxy; j++)
//         for (i = mregion -> minx; i <= mregion -> maxx; i++)
//           if (mregion -> included(i,j))
//           { //p = mregion->get_point(i,j);
//              //dis = mmodel->signed_distance(p);
//              //if (dis < 0) glColor3f(0.0,1.0,0.0);
//              //   else glColor3f(1.0,0.0,0.0);
//              //glVertex4f(p.x,p.y,p.z,mregion->theImage->normGL());
//              }    
//       glEnd();
//       glFlush();
//      }
//   }  
  
void description::print()
{ int i;
  mmodel->print();
  mregion->print();
  std::cout << "\nTotal error: " << total_error(*mregion,i) << "\n";
  }  
    
double& description::cameraX()
{  if (camera_set) return(cameraCoordinates[0]);
       else return(mregion->theImage->camOrientationX);
    }

double& description::cameraY()
{  if (camera_set) return(cameraCoordinates[1]);
       else return(mregion->theImage->camOrientationY);
    }

double& description::cameraZ()
{  if (camera_set) return(cameraCoordinates[2]);
       else return(mregion->theImage->camOrientationZ);
    }

void description::set_grow(int growing)
{ can_grow_var = growing; 
   }

int description::statOK(int *neigh, int *used, int nt, double thr, int &poz, int &neg)
{ int i,ok;

    for (i = poz+neg; i < nt; i++)
       if (used[neigh[i]])
      { if (mmodel->signed_distance(mregion->theImage->pixel(neigh[i])) < 0) neg++;
            else poz++;
         used[neigh[i]] = 0;
         }

   if (poz < neg) poz = neg;

   ok = oM.overflowP(nt,thr) > poz;

   return(ok);
   }

void description::analyze()
{ int i,k,ok,wr,poz,neg;
   int *neigh, *used, *accn, *accu;
   int nt, np = mregion->theImage->points(), ant;
   const int tries = 5;

  neigh = new int[np];
  used = new int[np];
  accn = new int[np];
  accu = new int[np];

   allowed->reset_all();
   memset(used,0,np<<2);
   nt = mregion->neighbourhood(neigh,used);


  for (i = 0; i < nt; i++)
     if (used[neigh[i]])
    { used[neigh[i]]=0;
       ant = 1;
       memset(accu,0,np<<2);
       accn[0] = neigh[i];
       accu[neigh[i]] = 1;
       poz = neg = 0;
       ant = mregion->neighbourhood(accn,accu, ant);
       ant = mregion->neighbourhood(accn,accu, ant);
       for (k = wr = 0; k < tries; k++) 
      { ant = mregion->neighbourhood(accn,accu, ant);
         if (!statOK(accn,accu,ant,0.99,poz,neg)) wr++;
         }
       
       ok = (wr <= (tries >> 1));
       if (ok) allowed->set_point(neigh[i]);
       }

   delete [] neigh;
   delete [] used;
   delete [] accn;
   delete [] accu;
   }

double description::residual2(region *r, int usual)
{ int i,j,k = discrepancy;
  double d, rez = 0.0;
  chiMSC *c1;

  if (usual) k = 0;
  if (r == NULL) r = mregion;
  switch (k)
  { case 0:
       rez = 0.0;
       for (i = r->minx; i <= r->maxx; i++)
          for (j = r->miny; j <= r->maxy; j++)
             if (r->included(i,j))
             { d = mmodel->signed_distance(r->get_point(i,j));
               rez += d*d;
               }
       break;
    case 1:
          c1 = new chiMSC(0.0,1.0,r->point_count());
          c1->fixate(this,r,0);
          d = c1->gprobability();
          if (d <= 0) rez = 1e10;
            else rez = -2.0*log(d);
          delete c1;
       break;
    }

  return(rez);
  }

double description::gprobability()
{ chiMSC *c1;
   double d,rez;

   c1 = new chiMSC(0.0,1.0,mregion->point_count());
   c1->fixate(this,mregion,0);
   d = c1->gprobability();
   if (d <= 0) rez = 1e10;
     else rez = -2.0*log(d);
  delete c1;
  return(rez);
  }

double description::MSC(description *d)
{ switch(msc)
  { case 0: return(MDL(d));
    case 1: return(AIC(d));
    case 2: return(BIC(d));
    }
  return(0.0);
  }

double description::MDL(description *d)
{ double rez,k1;
  
  k1 = 1.0;
  
  if (d == NULL || d == this)
  {
	std::cout << "MDL: " << mregion->point_count() << " " << residual2() << " " << mmodel->no_param() << "\n";
	  
	rez = k1 * mregion->point_count() - k2 * residual2() - k3 * mmodel->no_param();
    } else
    { region r(*mregion);
      r &= *d->mregion;
      rez = -k1*r.point_count()/2.0;
      }
  return(rez);
  }

double description::AIC(description *d)
{ double rez;

  if (d == NULL || d == this)
  { rez = 6.0 * mregion->point_count() - residual2() - 2 * mmodel->no_param();
    } else 
    { region r(*mregion);
      r &= *d->mregion;
      rez = - 3 * r.point_count();
      }
  return(rez);
  }

double description::BIC(description *d)
{ double rez;

  if (d == NULL || d == this)
  { rez = 3*log(2.0) * mregion->point_count() - residual2() - 
        log(1.0 * mregion->point_count()) * mmodel->no_param();
    } else
    { region r(*mregion);
      r &= *d->mregion;
      rez = -1.5*log(2.0)*r.point_count();
      }
  return(rez);
  }

description::~description()
{ delete mregion;
   delete mmodel;
   delete allowed;
  }

double description::stdev = 1.0;
int description::msc = 0;
int description::discrepancy = 0;
double description::k2 = 0.33;
double description::k3 = 0.33;

