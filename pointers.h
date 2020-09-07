#ifndef POINTERS_H
#define POINTERS_H

#include <cstdio>
#include <iostream>
#include <fstream>

#include "fftw3.h"
#include "confgel.h"

using namespace std;

class Pointers {

 public:

  Pointers(ConfGel *ptr) : 
    cgel(ptr),
    mem(ptr->mem),
    root(ptr->root), 
    rfall(ptr->rfall),
    confpg(ptr->confpg) {}

  virtual ~Pointers() {}

 protected:

  ConfGel *cgel;
  Memory *&mem;
  RootFns *&root;
  RootAll *&rfall;
  Confinpg *&confpg;

};



#endif
