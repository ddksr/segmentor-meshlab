// bitvect.C

#include "bitvect.h"
#include <stdio.h>


bitvect::bitvect(int n)
{ int i = 8*sizeof(unsigned int);
  log_bits_per_int = 0;
  while (i > 1)
  { i = i >> 1;
    log_bits_per_int++;
    }
  number = w = n;
  h = 1;
  bits = new unsigned int[size = 1 + ((n-1)>>log_bits_per_int)];
//  cout << "BCONS ";
//   for (i = 0; i < size; i++) bits[i] = 0;
  memset(bits,0,size<<2);
  }

bitvect::bitvect(int width, int height)
{ int i = 8*sizeof(unsigned int);
  log_bits_per_int = 0;
  while (i > 1)
  { i = i >> 1;
    log_bits_per_int++;
    }
  w = width;
  h = height;
  number = w * h;
  bits = new unsigned int[size = 1 + ((number-1)>>log_bits_per_int)];
//  cout << "BCONS ";
//  for (i = 0; i < size; i++) bits[i] = 0;
  memset(bits,0,size<<2);
  }

bitvect::bitvect(FILE *f)
{ int i = 8 * sizeof(unsigned int);
  log_bits_per_int = 0;
  while (i > 1)
  { i = i >> 1;
    log_bits_per_int++;
    }
  fread(&w,sizeof(int),1,f);
  fread(&h,sizeof(int),1,f);
  number = w * h;
  bits = new unsigned int[size = 1 + ((number-1)>>log_bits_per_int)];
//  cout << "BCONS ";
  fread(bits,sizeof(int),size,f);
  }  
  
bitvect::bitvect(bitvect& a)
{ log_bits_per_int = a.log_bits_per_int;
  w = a.w;
  h = a.h;
  number = a.number;
  bits = new unsigned int[size = a.size];
//  cout << "BCONS ";
//  for (i = 0; i < size; i++) bits[i] = a.bits[i];
  memcpy(bits,a.bits,size<<2);
  }

bitvect bitvect::operator=(const bitvect a)
{ 
  if (this!=&a)
  { delete [] bits;

    log_bits_per_int = a.log_bits_per_int;
    w = a.w;
    h = a.h;
    number = a.number;
    bits = new unsigned int[size = a.size];
//    for (i = 0; i < size; i++) bits[i] = a.bits[i];
    memcpy(bits,a.bits,size<<2);
    }
  return(*this);
  }

bitvect& bitvect::operator&=(const bitvect& a)
{ unsigned int i;
  if (w != a.w || h != a.h)
  { std::cout << "\nIncompatible types of bitvect in operator &=";
    } else
  { for (i = 0; i < size; i++) bits[i] &= a.bits[i];
    }
  return(*this);
  }

bitvect& bitvect::operator|=(const bitvect& a)
{ unsigned int i;
  if (w != a.w || h != a.h)
  { std::cout << "\nIncompatible types of bitvect in operator |=";
    } else
  { for (i = 0; i < size; i++) bits[i] |= a.bits[i];
    }
  return(*this);
  }

int bitvect::operator==(const bitvect& a)
{ unsigned int i;
  
  if (w != a.w || h != a.h) return(0);
  for (i = 0; i < size; i++)
    if (bits[i] != a.bits[i]) return(0);
  return(1);
  }


// returns 1 if a is subset of *this
int bitvect::subset(const bitvect& a)
{ unsigned int i;
  
  if (w != a.w || h != a.h)
  { std::cout << "\nIncompatible types of bitvect in function subset";
    return(0);
    } else
    { for (i = 0; i < number; i++)
        if (a.bit(i) && !bit(i)) return(0);
      }
  return(1);
  }

bitvect& bitvect::operator~()
{ unsigned int i;
  for (i = 0; i < size; i++) bits[i] = ~bits[i];
  return(*this);
  }

bitvect bitvect::operator&(const bitvect& a)
{ bitvect b = *this;
  b &= a;
  return(b);
  }

bitvect bitvect::operator|(const bitvect& a)
{ bitvect b = *this;
  b |= a;
  return(b);
  }

bitvect bitvect::dilation(int n, int *x, int *y)
{ unsigned int i,j;
  int k;
  unsigned int nx,ny;
  bitvect b(w,h);
  for (j = 0; j < h; j++)
    for (i = 0; i < w; i++)
      if (bit(i,j))
        for (k = 0; k < n; k++)
       	{ nx = i + x[k];
          ny = j + y[k];
          if (nx < w && ny < h)
            b.setbit(nx,ny);
	         }
  return(b);
  }

bitvect bitvect::erosion(int n, int *x, int *y)
{ unsigned int i,j;
  int k;
  unsigned int nx,ny;
  
  bitvect b(w,h);
 
  for (j = 0; j < h; j++)
    for (i = 0; i < w; i++)
    { b.setbit(i,j);
      for (k = 0; k < n && b.bit(i,j); k++)
      { nx = i + x[k];
        ny = j + y[k];
        if (nx < w && ny < h)
          if (!bit(nx,ny)) b.resetbit(i,j);
	}
      }
  return(b);
  }

  
std::ostream& operator<<(std::ostream& s, const bitvect& a)
{ int i,j,w = a.width(),h = a.height();
  for (j = 0; j < h; j++)
  { s << '\n';
    for (i = 0; i < w; i++)
      if (a.bit(i,j)) std::cout<<'1';
       	else std::cout<<'0';
    }
  return(s);
  }

void bitvect::save(FILE *f)
{ fwrite(&w,sizeof(int),1,f);
  fwrite(&h,sizeof(int),1,f);
  fwrite(bits,sizeof(int),size,f);
  }  
  
bitvect::~bitvect()
{ delete [] bits;
//  cout << "### ";
  }

/* How to use class bitvect

void main()
{ int i;
  bitvect b(5,5);

  b.setbit(0,0);
  b.setbit(1,0);
  b.setbit(4,1);
  cout << b;

  bitvect c = ~b;
  b.setbit(3,3);
  b.setbit(2,2);
  cout << c << (b & c) << (b | c);

  cin >> i;
  }
*/


