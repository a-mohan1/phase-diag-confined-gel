#include <complex>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "fftw3.h"
#include "confgel.h"
#include "confinpg.h"
#include "memory.h"
#include "rootfns.h"
#include "rootall.h"

using namespace std;

ConfGel::ConfGel() 
{
  mem = new Memory(this);
  root = new RootFns(this);
  rfall = new RootAll(this);
  confpg = new Confinpg(this);
}


ConfGel::~ConfGel() 
{ 
  delete confpg;
  delete rfall;
  delete root;
  delete mem;
}
