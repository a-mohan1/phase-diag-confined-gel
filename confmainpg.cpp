#include <complex>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "confinpg.h"
#include "confgel.h"
#include "fftw3.h"

using namespace std;

int main()
{
  ConfGel *confgel = new ConfGel();
  confgel->confpg->iterate();

  delete confgel;

  return 0;
}
