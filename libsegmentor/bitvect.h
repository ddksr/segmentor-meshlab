// bitvect.H

#ifndef LIBSEGMENTOR_BITVECT
#define LIBSEGMENTOR_BITVECT 1

#include <iostream>
#include <stdio.h>
#include <string.h>

class bitvect
{ public:
  unsigned int *bits;
  unsigned int w,h,size,log_bits_per_int;
  unsigned int number;

  bitvect(int n);
  bitvect(int width, int height);
  bitvect(FILE *f);
  bitvect(bitvect& a);                         // copy constructor
  bitvect operator=(const bitvect a);
  unsigned int bit(unsigned int i, unsigned int j) const
  { int z = i + j * w;
     int m = z >> log_bits_per_int;

    if (i < w && j < h)
    { if (bits[m] & (1<<(z-(m<<log_bits_per_int)))) return(1);
        else return(0);
      } else
    { std::cout << "\nbitvect element " << i << ',' << j << " out of range in bit";
      return(0);
      }
    }
    
  unsigned int bit(unsigned int i) const
  { int m = i >> log_bits_per_int;

    if (i < number)
    { return(bits[m] & (1<<(i-(m<<log_bits_per_int))));
      } else
    { std::cout << "\nbitvect element " << i << " out of range in bit";
      return(0);
      }
    }

  void setbit(unsigned int i, unsigned int j)
  { int z = i + j * w;
    int m = z >> log_bits_per_int;

  if (i < w && j < h) bits[m] |= (1 << (z - (m << log_bits_per_int)));
    else
     std::cout << "\nbitvect element " << i << ',' << j << " out of range in setbit";
    }
  void setbit(unsigned int i)
  { int m = i >> log_bits_per_int;

  if (i < number) bits[m] |= (1 << (i - (m << log_bits_per_int)));
    else
     std::cout << "\nbitvect element " << i << " out of range in setbit";
    }
  
  void resetbit(unsigned int i, unsigned int j)
  { int z = i + j * w;
    int m = z >> log_bits_per_int;

  if (i < w && j < h) bits[m] &= (~(1 << (z - (m << log_bits_per_int))));
    else
     std::cout << "\nbitvect element " << i << ',' << j << " out of range in setbit";
    }
  
  void resetbit(unsigned int i)
  { int m = i >> log_bits_per_int;

  if (i < number) bits[m] &= (~(1 << (i - (m << log_bits_per_int))));
    else
     std::cout << "\nbitvect element " << i << " out of range in setbit";
    }

  void reset()
  { memset(bits,0,size<<2);
    }

  void set()
  { int k = 0;
     k = ~k;
     memset(bits,k,size<<2);
     }
      
  bitvect& operator &=(const bitvect& a);
  bitvect& operator |=(const bitvect& a);
  int operator ==(const bitvect& a);
  int subset(const bitvect& a);
  unsigned int width() const { return(w); }
  unsigned int height() const { return(h); }
  bitvect& operator~();
  bitvect operator&(const bitvect& a);
  bitvect operator|(const bitvect& a);
  bitvect dilation(int n, int *x, int *y);
  bitvect erosion(int n, int *x, int *y);
  void save(FILE *f);
  ~bitvect();

  };

std::ostream& operator<<(std::ostream& s,const bitvect& a);

#endif

