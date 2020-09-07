#ifndef ROOTFNS_H
#define ROOTFNS_H

#include <fstream>
#include <cstdio>
#include <iostream>

#include "fftw3.h"
#include "pointers.h"

using namespace std;

class RootFns : protected Pointers {

 public:

  RootFns(class ConfGel *);
  ~RootFns() {};
  void lubksb(double**, int, int*, double*);
  void ludcmp(double**, int, int*, double&);
  void mnewt(int, double*, int, double, double, double);
  void fdjac(int, double*, double*, double**, double);
  void usrfun(double*, int, double*, double**, double);
  void lnsrch(int, double*, double, double*, double*, double*,
	      double&, double, bool&, double*, double);
  double MAX(double, double);
  void newt(double*, int, bool&, double);
  double fmin(double*, double*, int, double);
  void qrupdt(double**, double**, int, double*, double*);
  void qrdcmp(double**, int, double*, double*, bool&);
  void rotate(double**, double**, int, int, double, double);
  void rsolv(double**, int, double*, double*);
  double SIGN(double, double);
  void broydn(double*, int, bool&, double);

};



#endif
