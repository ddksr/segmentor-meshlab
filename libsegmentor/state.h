#ifndef LIBSEGMENTOR_STATE
#define LIBSEGMENTOR_STATE

struct LookupState {
  bool isSet;
  bool on;
  
  double a, b, c, e1, e2, kx, ky;
};

class State {
 public:
  static LookupState* lookup;
};

#endif
