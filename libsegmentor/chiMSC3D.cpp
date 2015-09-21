#include "chiMSC3D.h"

#include <math.h>
#include <stdio.h>

extern "C" {
#include "cdflib.h"
}

chiMSC::chiMSC(double m, double dev, int np, double thr)
{ int i,which,status;
   double dp = 0.0,bound,dq,mean,ddiv;
   int ndiv;

   ndiv = 5;
   if (ndiv < 5) ndiv = 5;
   ddiv = 1.0*ndiv;

   anyway = 0;
   which = 2;
   mean = m;
   stdev = dev;
   n = np;
   nmg = n / ndiv;
   conf = thr;
   param = 1;       
   if (nmg)
   { mg = new double[nmg];
     dist = new double[n];
     prob = new double[nmg];
     undodist = new double[n];
     cls = new int[nmg];
     undocls = new int[nmg];
     ndist = new int[n];
     undondist = new int[n];

     for (i = 0; i < n; i++) dist[i] = undodist[i] = 0.0;
     for (i = 0; i < nmg; i++) cls[i] = 0;
     for (i = 0; i < n; i++) ndist[i] = undondist[i] = 0;

     for (i = 1; i < nmg; i++)
     {  prob[i-1] = ddiv/n;
        dp += (ddiv/n);
        dq = 1.0-dp;
        cdfnor(&which,&dp,&dq,mg+i,&mean,&stdev,&status,&bound);
        //printf("\n%lf %lf",dp,mg[i]);
        }
     prob[nmg-1]=1-(nmg-1)*ddiv/n; 

     for (i = 0; i < n; i++) cls[classify(i)]++;
     for (i = 0; i < nmg; i++) 
     { undocls[i] = cls[i];
       //printf("\nDistribution:%lf %lf %d ",mg[i],prob[i],cls[i]);
       }
     }
   }

chiMSC::~chiMSC()
{  if (nmg)
   { delete [] mg;
     delete [] dist;
     delete [] undodist;
     delete [] cls;
     delete [] undocls;
     delete [] prob;
     delete [] ndist;
     delete [] undondist;
     }
   }

chiMSC::chiMSC(const chiMSC &c)
{ int i;
  anyway = c.anyway;
  stdev = c.stdev;
  mean = c.mean;
  conf = c.conf;
  n = c.n;
  nmg = c.nmg;
  param = c.param;
  if (nmg)
  { mg = new double[nmg];
    dist = new double[n];
    prob = new double[nmg];
    undodist = new double[n];
    cls = new int[nmg];
    undocls = new int[nmg];
    ndist = new int[n];
    undondist = new int[n];
    for (i = 0; i < nmg; i++)
    { mg[i] = c.mg[i];
      prob[i] = c.prob[i];
      cls[i] = c.cls[i];
      undocls[i] = c.undocls[i];
      }
    for (i = 0; i < n; i++)
    { dist[i] = c.dist[i];
      undodist[i] = c.undodist[i];
      ndist[i] = c.ndist[i];
      undondist[i] = c.undondist[i];
      }
    }
  }

double chiMSC::calcChi()
{  double chi = 0.0;
   int i;
   for (i = 0; i < nmg; i++) chi += (cls[i]-n*prob[i])*(cls[i]-n*prob[i]) / (n * prob[i]);
   return(chi);
   }

int chiMSC::classify(int i)
{ int l = 1,r = nmg-1,ni;

   if (nmg)
   { if (dist[i] < mg[l]) return(0);
        else if (dist[i] > mg[r]) return(nmg-1);
           else while (r-l > 1)
           {  ni = (l+r)/2;
              if (dist[i] > mg[ni]) l = ni;
                else if (dist[i] < mg[ni]) r = ni;
                  else return(ni);
              }
     }
   return(l);
   }

double chiMSC::absEvaluation(double &thr, int &r)
{ int which,status;
  double df,p,q,bound;
  if (nmg)
  { r = nmg - 1 - param;
    which = 2;
    p = conf;
    q = 1-conf;
    df = 1.0 * r;
    cdfchi(&which,&p,&q,&thr,&df,&status,&bound);
    return(calcChi());
    } else return(0.0);
  }

double chiMSC::probability()
{ double x = calcChi();
  int which,status;
  double p,q,df,bound,d;

  df = 1.0 * (nmg - 1 - param);
  which = 1;
  cdfchi(&which,&q,&p,&x,&df,&status,&bound);
  //printf("\nCLS: ");
  //for (i = 0; i < nmg; i++) printf("(%d-%lf) ",cls[i],n*prob[i]);
  
  if (status) d = -1e100;
    else if (p > 0.0) d = p;
       else if (df < x) d = df-x;
          else d = -1e100;
//  printf("%lf",d);
//  fflush(stdout);
  return(d);
  }

double chiMSC::gprobability()
{ double x = calcChi();
  int which,status;
  double p,q,df,bound,d,xf;

  xf = 1.0 * (nmg - 1 - param);
  if (xf < 1.0) df = 1.0;
    else df = xf;
  which = 1;
  cdfchi(&which,&q,&p,&x,&df,&status,&bound);
  //printf("\nCLS: ");
  //for (i = 0; i < nmg; i++) printf("(%d-%lf) ",cls[i],n*prob[i]);
  
  if (status) d = 0.0;
    else d = p * exp(xf-df);    // to penalize the leak of degrees of freedom
     
//  printf("%lf",d);
//  fflush(stdout);
  return(d);
  }

void chiMSC::clear()
{ int i,k;

  if (nmg)
  { for (i = 0; i < n; i++) dist[i] = undodist[i] = 0.0;
    for (i = 0; i < nmg; i++) cls[i] = undocls[i] = 0;
    for (i = 0; i < n; i++)
    { k = classify(i);
      cls[k]++;
      undocls[k]++;
      ndist[i] = undondist[i] = 0;
      }
    param = 1;
    anyway = 0;
    }
  }

void chiMSC::include(description *m, region *r, int fixate, int offset)
{ int i,j,k1,k2,k = 0;
  double d;

  if (nmg)
  { param = param+m->mmodel->no_param()+1;
    for (j = r->miny; j <= r->maxy; j++)
      for (i = r->minx; i <= r->maxx; i++)
        if (r->included(i,j))
        { d = m->mmodel->signed_distance(r->get_point(i,j));                  
          k1 = classify(k-offset);
          cls[k1]--;
          if (!ndist[k-offset]) dist[k-offset] = d;
            else dist[k-offset] = 0.0;
          ndist[k-offset]++;
          k2 = classify(k-offset);
          cls[k2]++;
          if (fixate)
          { undodist[k-offset] = dist[k-offset];
            undocls[k1]--;
            undocls[k2]++;
            undondist[k-offset]++;
            }
          k++;
          }
     }
  }

void chiMSC::exclude(description *m, region *r, int offset)
{ int i,j,k = 0;

  if (nmg)
  { for (j = r->miny; j <= r->maxy; j++)
      for (i = r->minx; i <= r->maxx; i++) 
      { dist[k-offset] = undodist[k-offset];
        ndist[k-offset] = undondist[k-offset];
        k++;
        }
    for (i = 0; i < nmg; i++) cls[i] = undocls[i];
    
    param = param - 1 - m->no_param();
    }
  }

void chiMSC::fixate(description *m, region *r, int offset)
{ if (nmg) include(m,r,1,offset);
  }


