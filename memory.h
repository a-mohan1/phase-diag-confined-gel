#ifndef MEMORY_H
#define MEMORY_H

#include <cstdio>
#include <iostream>
#include <fstream>

#include "fftw3.h"
#include "pointers.h"

using namespace std;

class Memory : protected Pointers {

 public:

  Memory(class ConfGel *);
  ~Memory() {};

  double **create_2d_double_array(int, int);
  void destroy_2d_double_array(double **);

};



#endif
