// linegeo.hxx


#ifndef _LINEGEO_H
#define _LINEGEO_H


#include "sline.h"

class h_sline;

class sh_sline;

class h_sline_set{
public:

  h_sline_set(int n);
  ~h_sline_set();
  const sh_sline * get_set(int j, int h) const;
  void put_into_set(int j, int h, sh_sline *hs);
  int get_k() const;
  void reset(int nn=0);

private:
  int sl_k;
  sh_sline   ***hls;

};



class h_sline : public sline{
  // An element in an ordered set of lines
public:
  h_sline(const sline &sl);
  ~h_sline();
  void put_before(h_sline *bef);
  int average( sline *asl,
	       int *nn,
	       double *sdevdir,
	       double *sdevnp) const;
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

  h_sline * get_next() const;
  void set_next( h_sline *);

  h_sline * get_prev() const;
  void set_prev( h_sline *);

private:
  h_sline *next;
  h_sline *prev;
};



class sh_sline : public h_sline{
  // An element in an ordered set of lines which is put in a mesh
  // ordered by direction
public:
  static h_sline_set hset;
  static void reset_hset(int nn=0);

  sh_sline(const h_sline &hsl);
  ~sh_sline();
  static  const sh_sline * init_cell_sequence();
  static  const sh_sline *get_next_filled_cell();
private:
  static int    shli;  
  static int    shlj;  
  void put_to_set_by_geo();

};


typedef sh_sline *  Ptr_sh_sline;
typedef sh_sline ** PtrPtr_sh_sline;


class oh_sline : public h_sline{
public:
  oh_sline(const h_sline &hsl,
	   int icn  =  1,
	   double idn= 0.0,
	   double idd = 0.0);
  ~oh_sline();
  void insert(oh_sline **beg);

  int         get_cn() const;
  double      get_dn() const;
  double      get_dd() const;

private:
  int    cn;
  double dn;
  double dd;
};




inline int h_sline_set::get_k() const
{
  return sl_k;
}

inline int        oh_sline::get_cn() const
{
  return cn;
}

inline double      oh_sline::get_dn() const
{
  return dn;
}

inline double      oh_sline::get_dd() const
{
  return dd;
}

#endif
