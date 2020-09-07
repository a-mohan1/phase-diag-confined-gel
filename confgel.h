#ifndef CONFGEL_H
#define CONFGEL_H

#include <cstdio>
#include <fstream>
#include <iostream>

#include "fftw3.h"

using namespace std;

class ConfGel {

 public:

  class Memory *mem;
  class Confinpg *confpg;
  class RootFns *root;
  class RootAll *rfall;

  ConfGel();
  ~ConfGel();

};


#endif
