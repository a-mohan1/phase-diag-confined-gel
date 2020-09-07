#ifndef CONFINPG_H
#define CONFINPG_H

#include <cstdio>
#include <iostream>
#include <fstream>

#include "fftw3.h"
#include "pointers.h"

using namespace std;

static const double PI = 3.14159265358979323846; 

class Confinpg : protected Pointers {

 public:

  Confinpg(class ConfGel *);
  ~Confinpg();
  void readinp();
  void iterate();
  void iterF();
  double SolveF1(double**, double**, int, double, int, int);
  void calcphi();
  void calcphisol();
  double calcavphi();
  void diffsolve_0();
  void diffsolve_1s(int);
  void diffsolve_1t(int);
  void diffsolve_1t0();
  void diffeqsolver(int, double*, double*);
  void calcalpha();
  void calcalphasol();
  void calcmoments();
  void calcmomentsol();
  void calcphiumon();
  void calcF0s();
  void calcF0t();
  void printres();
  void calcfns(double*, double*, double);
  void calcallfns(double*, double*);

 private:  

  double len;
  double uparam;
  double h;
  double NA;
  double f;

  double LamF1;
  double rootfac;

  int Nx;
  int Ns;
  int Ntheta;

  double zA;

  int wAintype;
  int Fsscheme, Ftscheme;

  double hparam;
  double DelTheta;
  double Dels;

  int MAXITER;
  double MAXERRF;

  int numarq0, numarq1;

  double Nn, Nw, Nnsol, Nwsol;

  double *inout;
  fftw_complex *infwd, *outfwd;
  fftw_complex *inbkd, *outbkd;
  fftw_plan fftwp;
  fftw_plan fftwdftfp, fftwdftbp;

  double **F1s, **F1t, **F0s, **F0t;
  double **q1s, **q1t;
  double *q0;
  double *q1t0;
  double *wA;
  double *phitot, *phisol, *phie, *phiue, *phiesol, *phiuesol;
  double *alpha, *alphasol, *phiunrmon;

  double *fyn1, *fyn2;

  double hval;
};



#endif
