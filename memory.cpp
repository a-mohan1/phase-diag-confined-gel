#include <complex>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "fftw3.h"
#include "memory.h"

using namespace std;

Memory::Memory(ConfGel *cgel) : Pointers(cgel) {}

double **Memory::create_2d_double_array(int n1, int n2)
{
  double *data = (double *) malloc(n1*n2*sizeof(double));
  double **array = (double **) malloc(n1*sizeof(double *));

  int n = 0;
  for (int i = 0; i < n1; i++) {
    array[i] = &data[n];
    n += n2;
  }

  return array;
}


void Memory::destroy_2d_double_array(double **array)
{
  if (array == NULL) return;
  free(array[0]);
  free(array);
}
