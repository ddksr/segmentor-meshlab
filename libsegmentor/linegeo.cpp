// linegeo.cxx
// 9-May-97
// GL Cardiff
//

#ifndef LIBSEGMENTOR_LINEGEO
#define LIBSEGMENTOR_LINEGEO

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include <iostream>
#include "rcad_vector2.h"
#include "linegeo.h"

#define M_PI 3.14159265358979323846

h_sline_set sh_sline::hset(8);
int         sh_sline::shli(-1);
int         sh_sline::shlj(-1);

h_sline_set::h_sline_set(int k):
sl_k(0),
hls(NULL)
{
  // Produces an empty array of k*(k+1)/2 pointers

  int i,j;

  if (k > 0) {
    hls = new PtrPtr_sh_sline[k];
    sl_k = k;

    for (i=0; i < k; ++i) {
      hls[i] = new Ptr_sh_sline[i+1];
      for (j=0; j <= i; ++j) {
	hls[i][j] = NULL;
      }
    }
  }
}


h_sline_set::~h_sline_set()
{

  int i,j;

  for (i=0; i < sl_k; ++i) {
    for(j=0; j <= i; ++j) {
      delete hls[i][j];
    }
    delete [] (hls[i]);
  }

  delete [] hls;

  sl_k=0;
  hls = NULL;

}




const sh_sline *
h_sline_set::get_set(int j, int h) const
{
  // Returns the header pointer of the set

  if ((j < sl_k) && (h <= j) && (h >= 0))
    return hls[j][h];
  else
    return (sh_sline *) (NULL);
}



void
h_sline_set::reset(int nn)
{
  int i,j;

  if (nn > 0) {
    for (i=0; i < sl_k; ++i) {
      for(j=0; j <= i; ++j) {
	delete hls[i][j];
      }
      delete [] (hls[i]);
    }

    delete [] hls;
  
    hls = new PtrPtr_sh_sline[nn];
    sl_k = nn;

    for (i=0; i < nn; ++i) {
      hls[i] = new Ptr_sh_sline[i+1];
      for (j=0; j <= i; ++j) {
	hls[i][j] = NULL;
      }
    }
  } else {
    for (i=0; i < sl_k; ++i) {
      for(j=0; j <= i; ++j) {
	delete hls[i][j];
	hls[i][j] = NULL;
      }
    }
  }
}




void
h_sline_set::put_into_set(int j, int h, sh_sline *hs)
{

  if ((j < sl_k) && (h <= j) && (h >= 0)) {
    if (hls[j][h] != NULL) {
      hls[j][h]->put_before(hs);
    }
    hls[j][h] = hs;
  }
}

// h_sline:


h_sline::h_sline(const sline &nsl):
sline(nsl),
next(NULL),
prev(NULL)
{
}

h_sline::~h_sline()
{
  delete next;
  next = NULL;
}

void
h_sline::put_before(h_sline *bef)
  // bef should not be connected!
  // Previously:   prev=>this=>next; null=>bef=>null
  // Subsequently:     prev=>bef=>this=>next
{
  h_sline *x;

  x = prev;
  prev = bef;
  bef->prev = x;
  if (x != NULL) x->next   = bef;
  bef->next = this;
}



int h_sline::average(
		      sline *asl,
		      int *nn,
		      double *sdevdir,
		      double *sdevnp) const
  // Computes the average separately for the finite and infinite solutions
  // finite sol  : asl[0], nn[0], sdevdir[0], *sdevnp
  // infinite sol: asl[1], nn[1], sdevdir[1]
  // asl is the line
  // nn  is the cardinality in the list
  // sdevdir is standard deviation of the direction
  // sdevnp is standard deviation of the nearest point to the origin
  // (this last one is not present at the infinite solutions).
  // Both 'sdevdir' and 'sdevnp' can be NULL meaning that the particular
  // deviation is not asked for.
  // The function returns the sum of the lines in the list.
{
  int    n;
  double a0dir[3];
  double a1dir[3];
  double anp[3];
  double dir[3];
  double np[3];
  const h_sline  *h;
  double dd[2],
         dn;

  for (n=0; n < 3; ++n)
    a0dir[n] = a1dir[n] = anp[n] = 0.0;

  nn[0] = nn[1] = 0;
  for (h = this; h != NULL; h = h->next) {
    h->direction(dir);
    if ((dir[2] < 0.0) ||
	((dir[2] == 0.0) && (dir[1] < 0.0)) ||
	((dir[2] == 0.0) && (dir[1] = 0.0) && (dir[0] < 0.0)))
      mult(dir,  (double) -1.0, dir);

    if (h->finite()) {
      h->nearp(np);
      add( a0dir, a0dir, dir);
      add( anp,  anp, np);
      ++nn[0];
    } else {
      add( a1dir, a1dir, dir);
      ++nn[1];
    }
  }

  if (nn[0] > 0) {
 //   mult(anp, 1.0/((double) n), anp);

// possible error. ADM thinks div should be nn[0]
mult(anp, 1.0/((double) nn[0]), anp);


    unit(a0dir);
    asl[0] = asl[1] = sline(anp, a0dir);
  }

  if (nn[1] > 0) {
    unit(a1dir);
    asl[1] = sline(NULL, a1dir);
    if (nn[0]==0) {
      asl[0] = asl[1];
    }
  }
  n = nn[0]+nn[1];

  if ((sdevdir != NULL) || (sdevnp != NULL)) {

    dd[0] = dd[1] = dn = 0.0;

    if (n > 1) {
      for (h = this; h != NULL; h = h->next) {
	h->direction(dir);
	if ((dir[2] < 0.0) ||
	    ((dir[2] == 0.0) && (dir[1] < 0.0)) ||
	    ((dir[2] == 0.0) && (dir[1] = 0.0) && (dir[0] < 0.0)))
	  mult(dir, -1.0, dir);

	if (h->finite()) {
	  sub( dir, a0dir, dir);
	  dd[0] += s_mult(dir,dir);
	  h->nearp(np);
	  sub( np,  anp, np);
	  dn += s_mult(np,np);
	} else {
	  sub( dir, a1dir, dir);
	  dd[1] += s_mult(dir,dir);
	}
      }
      if (nn[0] > 1) {
	dd[0] = sqrt ( dd[0] / ((double)(nn[0]-1)));
	dn = sqrt ( dn / ((double)(nn[0]-1)));
      }
      if (nn[1] > 1) {
	dd[1] = sqrt ( dd[1] / ((double)(nn[1]-1)));
      }
    }
    if (sdevdir != NULL) {
      sdevdir[0] = dd[0];
      sdevdir[1] = dd[1];
    }
    if (sdevnp != NULL) {
      *sdevnp = dn;
    }
  }

  return n;

}


h_sline * h_sline::get_next() const
{
  return next;
}

void h_sline::set_next( h_sline * nx)
{
  next = nx;
}

h_sline * h_sline::get_prev() const
{
  return prev;
}

void h_sline::set_prev( h_sline *pr)
{
  prev = pr;
}



// sh_sline:


sh_sline::sh_sline(const h_sline &hs):
h_sline(hs)
{
  // Nyu'zni kell = Should be newn

  put_to_set_by_geo();
}

sh_sline::~sh_sline()
{
}

void
sh_sline::reset_hset(int nn)
{
  hset.reset(nn);
  shli = shlj = -1;
}


void
sh_sline::put_to_set_by_geo()
{
  int k=hset.get_k();
  double dir[3];
  double phi;
  double c;
  int j,h;

  direction(dir);
  if (dir[2] < 0.0) {
    mult(dir, -1.0, dir);
    reverse();
  }

  if (dir[2] >= 1.0) {
    dir[2] = 1.0;
    phi = 0.0;
  } else {
    phi = atan2(dir[1],dir[0]);
  }

  if ((dir[2] == 0.0) && (phi >= M_PI)) {
    reverse();
    phi -= M_PI;
  }

  c = 4.0*(1.0-dir[2])*((double) (k*(k+1)));

  j = (int) (0.5 *(1.0 + sqrt(1.0+c)));

  if (j > k) j = k;

  c = 0.5*((double) j) * phi/M_PI;

  --j;
  h = (int) c;

  if (h > j) h = j;

  hset.put_into_set(j,h,this);
}


const sh_sline * sh_sline::init_cell_sequence()
{
  int k=hset.get_k();
  const sh_sline * q;

  shli = shlj = 0;
  q = hset.get_set(0,0);
  while (q == NULL) {
    if (shlj >= shli) {
      ++shli;
      shlj=0;
      if (shli >= k) {
	shli = 0;
	break;
      }
    } else {
      ++shlj;
    }
    q = hset.get_set(shli,shlj);
  }
  return q;
}
       

const sh_sline *sh_sline::get_next_filled_cell()
{
  int k=hset.get_k();
  const sh_sline * q = 0;

  if ((shli < 0) || (shlj < 0)) {
    q = init_cell_sequence();
  } else {
    do {
      if (shlj >= shli) {
	++shli;
	shlj=0;
	if (shli >= k) {
	  shli = 0;
	  break;
	}
      } else {
	++shlj;
      }
      q = hset.get_set(shli,shlj);
    } while (q == NULL);
  }

  return q;
}

// oh_sline


oh_sline::oh_sline(const h_sline &hsl,
	 int icn,
	 double idn,
	 double idd):
h_sline(hsl),
cn(icn),
dn(idn),
dd(idd)
{
}

oh_sline::~oh_sline()
{
}

void oh_sline::insert(oh_sline **beg)
{
  oh_sline *x;
  oh_sline *y;

  y = NULL;
  x = *beg;

  while (x != NULL) {
    if (((x)->get_cn() < cn) ||
	(((x)->get_cn() == cn) && 
	 (dn < (x)->get_dn()))) {
      x->put_before(this);
      if (x == *beg) {
	*beg = this;
      }
      break;
    }
    y = x;
    x = (oh_sline *) y->get_next();
  }

  if (x == NULL) {
    set_prev(y);
    if (y != NULL) {
      y->set_next(this);
    } else {
      *beg = this;
    }
  }
}



#endif
