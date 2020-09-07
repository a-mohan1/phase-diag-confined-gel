#include <complex>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "fftw3.h"
#include "confinpg.h"
#include "memory.h"
#include "rootfns.h"
#include "rootall.h"

using namespace std;

Confinpg::Confinpg(ConfGel *cgel) : Pointers(cgel)
{
  // read input file
  readinp();  

  numarq0 = Nx*(Ns+1);
  numarq1 = Nx*(2*Ns+1);

  // fftw DST-I plan
  inout = (double *) fftw_malloc(Nx*sizeof(double));
  fftwp = fftw_plan_r2r_1d(Nx, inout, inout, FFTW_RODFT00, FFTW_MEASURE);

  // fftw DFT plans
  infwd = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * Nx);
  outfwd = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * Nx);
  inbkd = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * Nx);
  outbkd = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * Nx);
  fftwdftfp = fftw_plan_dft_1d(Nx, infwd, outfwd, FFTW_FORWARD, FFTW_MEASURE);
  fftwdftbp = fftw_plan_dft_1d(Nx, inbkd, outbkd, FFTW_BACKWARD, FFTW_MEASURE);

  // allocate vars
  // smaller root (sol) and larger root (normalization) of F1
  F1s = mem->create_2d_double_array(Nx, Ntheta);
  F0s = mem->create_2d_double_array(Nx, Ntheta);
  q1s = mem->create_2d_double_array(numarq1, Ntheta-1);
  F1t = mem->create_2d_double_array(Nx, Ntheta);
  F0t = mem->create_2d_double_array(Nx, Ntheta);
  q1t = mem->create_2d_double_array(numarq1, Ntheta-1);
  q0 = (double *) malloc(numarq0*sizeof(double));
  q1t0 = (double *) malloc(numarq0*sizeof(double));
  wA = (double *) malloc(Nx*sizeof(double));  
  phitot = (double *) malloc(Nx*sizeof(double));
  phie = (double *) malloc(Nx*sizeof(double));
  phiue = (double *) malloc(Nx*sizeof(double));
  alpha = (double *) malloc(Nx*sizeof(double));
  phisol = (double *) malloc(Nx*sizeof(double));
  phiesol = (double *) malloc(Nx*sizeof(double));
  phiuesol = (double *) malloc(Nx*sizeof(double));
  alphasol = (double *) malloc(Nx*sizeof(double));
  phiunrmon  = (double *) malloc(Nx*sizeof(double));

  // Adams-Bashforth-Moulton vars
  fyn1 = (double *) malloc(Nx*sizeof(double));
  fyn2 = (double *) malloc(Nx*sizeof(double));

}


Confinpg::~Confinpg() 
{ 
  mem->destroy_2d_double_array(F1s);
  mem->destroy_2d_double_array(F0s);
  mem->destroy_2d_double_array(q1s);
  mem->destroy_2d_double_array(F1t);
  mem->destroy_2d_double_array(F0t);
  mem->destroy_2d_double_array(q1t);
  free(q0);
  free(q1t0);
  free(wA);
  free(phitot);
  free(phie);
  free(phiue);
  free(alpha);
  free(phisol);
  free(phiesol);
  free(phiuesol);
  free(alphasol);
  free(phiunrmon);

  free(fyn1);
  free(fyn2);

  // delete fftw plans
  fftw_destroy_plan(fftwp);
  fftw_destroy_plan(fftwdftfp);
  fftw_destroy_plan(fftwdftbp);
  fftw_free(inout);
  fftw_free(infwd);
  fftw_free(outfwd);  
  fftw_free(inbkd);
  fftw_free(outbkd);

}


void Confinpg::readinp()
{
  char tempstr[50];
  ifstream fin;
  fin.open("infile.dat", ios::in|ios::binary);

  while(fin >> tempstr) {
     if(strcmp(tempstr, "len")==0) {
       fin >> tempstr;
       len = atof(tempstr);
     }
     else if(strcmp(tempstr, "uparam")==0){
       fin >> tempstr;
       uparam = atof(tempstr);
     }
     else if(strcmp(tempstr, "h")==0){
       fin >> tempstr;
       h = atof(tempstr);
     }
     else if(strcmp(tempstr, "NA")==0){
       fin >> tempstr;
       NA = atof(tempstr);
     }
     else if(strcmp(tempstr, "f")==0){
       fin >> tempstr;
       f = atof(tempstr);
     }
     else if(strcmp(tempstr, "LamF1")==0){
       fin >> tempstr;
       LamF1 = atof(tempstr);
     }
     else if(strcmp(tempstr, "rootfac")==0){
       fin >> tempstr;
       rootfac = atof(tempstr);
     }
     else if(strcmp(tempstr, "Nx")==0){
       fin >> tempstr;
       Nx = atoi(tempstr);
     }
     else if(strcmp(tempstr, "Ns")==0){
       fin >> tempstr;
       Ns = atoi(tempstr);
     }
     else if(strcmp(tempstr, "Ntheta")==0){
       fin >> tempstr;
       Ntheta = atoi(tempstr);
     }
     else if(strcmp(tempstr, "zAinitial")==0){
       fin >> tempstr;
       zA = atof(tempstr);
     }
     else if(strcmp(tempstr, "wAinitial")==0){
       fin >> tempstr;
       wAintype = atoi(tempstr);
     }
     else if(strcmp(tempstr, "Fsscheme")==0){
       fin >> tempstr;
       Fsscheme = atoi(tempstr);
     }
     else if(strcmp(tempstr, "Ftscheme")==0){
       fin >> tempstr;
       Ftscheme = atoi(tempstr);
     }
     else if(strcmp(tempstr, "MAXITER")==0){
       fin >> tempstr;
       MAXITER = atoi(tempstr);
     }
     else if(strcmp(tempstr, "MAXERRF")==0){
       fin >> tempstr;
       MAXERRF = atof(tempstr);
     }
     else {
       cout<<"undefined input"<<endl;
       exit(1);
     }
  }

  fin.close();

  hparam = exp(h)/NA;
  DelTheta = 1.0/(double(Ntheta)-1.0);
  Dels = 1.0/double(Ns);
}


void Confinpg::iterate()
{

  // initialize field
  if (wAintype == 3) {
    // read input file
    int i = -1;
    char tempstr[50];
    ifstream wfin;
    wfin.open("wAin.dat", ios::binary|ios::in);
    while(wfin >> tempstr) {
       i++;
       wA[i] = atof(tempstr);
    }
    wfin.close();
    if(i != Nx-1) {
      cout<<"Error in wAin file"<<endl;
      exit(1);
    }
  }
  else {
    for(int i = 0; i < Nx; i++) {
      if (wAintype == 0) { 
	// homogeneous
	wA[i] = uparam;
      }
      else if (wAintype == 1) {
	// sine_av
	wA[i] = uparam*PI/2.0*sin(PI* double(i+1)/double(Nx+1));
      }
      else if (wAintype == 2) {
	// sine_max
	wA[i] = uparam*sin(PI* double(i+1)/double(Nx+1));
      }
      else {
	cout<<"unknown wAinitial"<<endl;
	exit(1);
      }
    }
  }

  // calculate F1, zA and wA
  iterF();
 
  // calculate and print results
  printres();

}


void Confinpg::iterF()
{
  int iterF;
  double errF;
  int indexs1, indexs2;    
  double *xvec = (double *) malloc((2*Nx+1)*sizeof(double));
  double *F1new = (double *) malloc(Nx*sizeof(double));
  double tolx, tolf;
  tolx = tolf = MAXERRF;
  int ntrial = MAXITER;  
  bool check = true;
  double theta;

  // initial guess
  diffsolve_0();
  for(int ix = 0; ix < Nx; ix++) {
    indexs1 = ix*(Ns+1) + Ns;
    F1t[ix][0] = pow(q0[indexs1], f-1.0);
    F1t[ix][Ntheta-1] = F1t[ix][0]*rootfac;
    xvec[ix] = F1t[ix][Ntheta-1];
    xvec[Nx+ix] = wA[ix];
  }
  xvec[2*Nx] = zA;

  // Root Finding
  // Newton-Raphson with line searches and back tracking
  rfall->newt(xvec, 2*Nx+1, check);
  if(check) {
    cout<<"root-finding failed"<<endl;
    exit(1);
  }

  for(int ix = 0; ix < Nx; ix++) {
    F1t[ix][Ntheta-1] = xvec[ix];
    wA[ix] = xvec[Nx+ix];
  }
  zA = xvec[2*Nx];
  diffsolve_0();
  diffsolve_1t(Ntheta-1);
  calcphi();
  free(xvec);

  // theta=0
  for(int ix = 0; ix < Nx; ix++) {
    indexs1 = ix*(Ns+1) + Ns;
    // theta=0:
    F1s[ix][0] = pow(q0[indexs1], f-1.0);
    F1t[ix][0] = pow(q0[indexs1], f-1.0);
    // initial condition and initial guess for theta=1, theta=DelTheta
    F1s[ix][Ntheta-1] = F1s[ix][0];
    F1s[ix][1] = F1s[ix][0];
  }

  // sol, theta=1
  errF = 1.0e6;
  iterF = 0;
  while(errF > MAXERRF) {
    // solve FK for theta=1 using theta=0
    // or previous iteration as initial condition
       
    diffsolve_1s(Ntheta-1);
    // this gives new q1s at theta=1

    // update F1s
    errF = SolveF1(F1s, q1s, Ntheta-1, LamF1, iterF, Fsscheme);

    // increment iteration counter
    iterF++;
    if(iterF > MAXITER) {
      cout<<"Max F iterations exceeded"<<endl;
      exit(1);
    }
    // iterate until convergence
  }
  // end of F1s (x, theta=1) calculation

  // solve for F1s,t for remaining theta in (0,1)
  for(int ithe = 1; ithe < Ntheta-1; ithe++) {

    // solve for F1s for each theta>0 until convergence
    errF = 1.0e6;
    iterF = 0;
    while(errF > MAXERRF) {
      // solve FK using previous theta 
      // or previous iteration as initial condition
       
      diffsolve_1s(ithe);
      // this gives new q1s at current theta

      // update F1s
      errF = SolveF1(F1s, q1s, ithe, LamF1, iterF, Fsscheme);

      // increment iteration counter
      iterF++;
      if(iterF > MAXITER) {
	cout<<"Max F iterations exceeded"<<endl;
	exit(1);
      }
      // iterate until convergence
    }

    // initial guess for next theta
    if(ithe < Ntheta-2) {
      for(int ix = 0; ix < Nx; ix++) {
	F1s[ix][ithe+1] = F1s[ix][ithe];
      }
    }

    // solve for F1t for each theta>0 until convergence
    // q1 is calculated during root-finding 
   
    // initial guess
    for(int ix = 0; ix < Nx; ix++) {
	F1t[ix][ithe] = F1s[ix][ithe] * rootfac;
    }

    theta = double(ithe)*DelTheta;

    // initial guess
    for(int ix = 0; ix < Nx; ix++) {
      F1new[ix] = F1t[ix][ithe];
    }

    // Root Finding
    if(Ftscheme == 0) {
       // Newton-Raphson
       root->mnewt(ntrial, F1new, Nx, tolx, tolf, theta);
    } else if (Ftscheme == 1) {
       // Newton-Raphson with line searches and back tracking
       root->newt(F1new, Nx, check, theta);
       if(check) {
	 cout<<"NR with line searches and backtracking failed"<<endl;
	 exit(1);
       }
    } else if (Ftscheme == 2) {
       // Broyden's method
       root->broydn(F1new, Nx, check, theta);
       if(check) {
	 cout<<"Broyden's method failed"<<endl;
	 exit(1);
       }
    } else {
       cout<<"unknown F1tscheme"<<endl;
       exit(1);
    }

    // update F1t at theta
    for(int ix = 0; ix < Nx; ix++) {
      F1t[ix][ithe] = F1new[ix];
    }

    diffsolve_1t(ithe);
    // this gives new q1t at theta

    // increment theta
  }
  free(F1new);


  // write F1s,t for all theta to file
  ofstream fout;
  fout.open("F1theta.dat", ios::binary|ios::out);
  for(int ithe = 0; ithe < Ntheta; ithe++) {
    fout<<ithe<<endl;
    for(int ix = 0; ix < Nx; ix++) {
      fout<<ix<<"  "<<F1s[ix][ithe]<<"  "<<F1t[ix][ithe]<<endl;
    }
  }
  fout.close();

  // write error to file
  double werr, F1t1err;
  fout.open("errfile.dat", ios::binary|ios::out);
  double avphi = calcavphi();
  fout<<1.0-avphi<<endl;
  for(int ix = 0; ix < Nx; ix++) {
      // q0 (s=N)
      indexs1 = ix*(Ns+1) + Ns;
      // q1 (s=2N)
      indexs2 = ix*(2*Ns+1) + 2*Ns;
      F1t1err = pow(q0[indexs1] + 
              zA*hparam*f*q1t[indexs2][Ntheta-2], f-1.0) 
               - F1t[ix][Ntheta-1];
      werr = -wA[ix]/uparam + phitot[ix];
      fout<<ix<<"  "<<F1t1err<<"  "<<werr<<endl;
  }
  fout.close();


}


double Confinpg::SolveF1(double **F1mat, double **q1ar, 
			 int ithe, double LamF, 
                         int iternum, int Fscheme) 
{
  // update F1s,t for theta at ithe
 
  double theta = double(ithe)*DelTheta;
  double errorF = 0.0;
  int indexs1, indexs2;

  if(Fscheme == 2) {
    // semi-implicit 
    // with linear term added and subtracted
    double km;
    fftw_complex *F1rold, *term2r;
    fftw_complex *F1kold, *term2k;
    fftw_complex *F1rnew, *F1knew;
    F1rold = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*Nx);
    term2r = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*Nx);
    F1kold = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*Nx);
    term2k = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*Nx);  
    F1rnew = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*Nx);
    F1knew = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*Nx);

    for(int ix = 0; ix < Nx; ix++) {
      // q0 (s=N)
      indexs1 = ix*(Ns+1) + Ns;
      // q1 (s=2N)
      // q1: theta_index=0:ithe=1
      indexs2 = ix*(2*Ns+1) + 2*Ns;
      // Re
      F1rold[ix][0] = F1mat[ix][ithe];
      term2r[ix][0] = pow(q0[indexs1] + 
            zA*hparam*f*theta*q1ar[indexs2][ithe-1], f-1.0);
      // Im
      F1rold[ix][1] = 0.0;
      term2r[ix][1] = 0.0;
    }

    fftw_execute_dft(fftwdftfp, F1rold, F1kold);
    fftw_execute_dft(fftwdftfp, term2r, term2k);

    for(int m = 0; m < Nx; m++) {
      if(m <= Nx/2) {
	km = 2.0*PI*double(m)/len;
      } else {
	km = 2.0*PI*double(m-Nx)/len;
      }

      for(int ir = 0; ir < 2; ir++) {

        F1knew[m][ir] = F1kold[m][ir] + LamF*term2k[m][ir] 
        - LamF*(f-1.0)*exp(-uparam*f)*zA*hparam*f*theta*exp(-2.0*km*km)
            *F1kold[m][ir];
        F1knew[m][ir] /= (1.0 + LamF 
         - LamF*(f-1.0)*exp(-uparam*f)*zA*hparam*f*theta*exp(-2.0*km*km));
      }
    }

    fftw_execute_dft(fftwdftbp, F1knew, F1rnew);
 
    // update F1mat[ix][ithe] and calculate error
    // normalization by N_x
    for(int ix = 0; ix < Nx; ix++) {
      F1rnew[ix][0] /= double(Nx);
      // Re
      // discard Im
      errorF += fabs(F1rnew[ix][0] - F1mat[ix][ithe]);
      F1mat[ix][ithe] = F1rnew[ix][0];
    }

    fftw_free(F1rold);
    fftw_free(term2r);
    fftw_free(F1kold);
    fftw_free(term2k);
    fftw_free(F1rnew);
    fftw_free(F1knew);

  }
  else if (Fscheme == 3) {

    // Adams-Bashforth-Moulton scheme

    double F1new;
    double funcval, fypred;
    int Ncs = 2*Ns+1;
    double *ypred = (double *) malloc(Nx*sizeof(double));
    double *qpred = (double *) malloc(Nx*Ncs*sizeof(double));

    if (iternum == 0) {
      for(int ix = 0; ix < Nx; ix++) {
        // q0 (s=N)
        indexs1 = ix*(Ns+1) + Ns;
        // q1 (s=2N)
        // q1: theta_index=0:ithe=1
        indexs2 = ix*(2*Ns+1) + 2*Ns;

        funcval = pow(q0[indexs1] + 
	zA*hparam*f*theta*q1ar[indexs2][ithe-1], f-1.0) 
          - F1mat[ix][ithe];

        F1new = F1mat[ix][ithe] + LamF*funcval;
        fyn1[ix] = funcval;

        // compute error in F1
        errorF += fabs(F1new - F1mat[ix][ithe]);

	F1mat[ix][ithe] = F1new;

        cout<<pow(q0[indexs1] +
              zA*hparam*f*theta*q1ar[indexs2][ithe-1], f-1.0)
              -  F1mat[ix][ithe]<<"  ";

      }
    } else if (iternum == 1) {
      for(int ix = 0; ix < Nx; ix++) {
        // q0 (s=N)
        indexs1 = ix*(Ns+1) + Ns;
        // q1 (s=2N)
        // q1: theta_index=0:ithe=1
        indexs2 = ix*(2*Ns+1) + 2*Ns;
        funcval = pow(q0[indexs1] + 
	zA*hparam*f*theta*q1ar[indexs2][ithe-1], f-1.0) 
          - F1mat[ix][ithe];

	fyn2[ix] = fyn1[ix];
	fyn1[ix] = funcval;
	F1new = F1mat[ix][ithe] + LamF*funcval;

        // compute error in F1
        errorF += fabs(F1new - F1mat[ix][ithe]);

	F1mat[ix][ithe] = F1new;

        cout<<pow(q0[indexs1] +
              zA*hparam*f*theta*q1ar[indexs2][ithe-1], f-1.0)
              -  F1mat[ix][ithe]<<"  ";
      }
    } else {
      for(int ix = 0; ix < Nx; ix++) {
        // q0 (s=N)
        indexs1 = ix*(Ns+1) + Ns;
        // q1 (s=2N)
        // q1: theta_index=0:ithe=1
        indexs2 = ix*(2*Ns+1) + 2*Ns;

        funcval = pow(q0[indexs1] + 
	zA*hparam*f*theta*q1ar[indexs2][ithe-1], f-1.0) 
          - F1mat[ix][ithe];

	ypred[ix] = F1mat[ix][ithe] + LamF/12.0*
	  (23.0*funcval - 16.0*fyn1[ix] + 5.0*fyn2[ix]);
      }

      diffeqsolver(Ncs, ypred, qpred);	

      for(int ix = 0; ix < Nx; ix++) {
        // q0 (s=N)
        indexs1 = ix*(Ns+1) + Ns;
        // q1 (s=2N)
        // q1: theta_index=0:ithe=1
        indexs2 = ix*(2*Ns+1) + 2*Ns;

        funcval = pow(q0[indexs1] + 
	zA*hparam*f*theta*q1ar[indexs2][ithe-1], f-1.0) 
          - F1mat[ix][ithe];

        fypred = pow(q0[indexs1] + 
	zA*hparam*f*theta*qpred[indexs2], f-1.0) - ypred[ix];

        F1new = F1mat[ix][ithe] + LamF/12.0*
	  (5.0*fypred + 8.0*funcval - fyn1[ix]);

        fyn2[ix] = fyn1[ix];
	fyn1[ix] = funcval;

        // compute error in F1
        errorF += fabs(F1new - F1mat[ix][ithe]);

	F1mat[ix][ithe] = F1new;

        cout<<pow(q0[indexs1] +
              zA*hparam*f*theta*q1ar[indexs2][ithe-1], f-1.0)
              -  F1mat[ix][ithe]<<"  ";


      }

    }

    cout<<endl;

    free(ypred);
    free(qpred);

  
  } 
  else if (Fscheme == 4) {

    // second order Runge-Kutta method

    double F1new;
    double funcval;
    int Ncs = 2*Ns+1;
    double *yint, *qyint;
    yint = (double *) malloc(Nx*sizeof(double));
    qyint = (double *) malloc(Nx*Ncs*sizeof(double));
  
    for(int ix = 0; ix < Nx; ix++) {
      // q0 (s=N)
      indexs1 = ix*(Ns+1) + Ns;
      // q1 (s=2N)
      // q1: theta_index=0:ithe=1
      indexs2 = ix*(2*Ns+1) + 2*Ns;
      funcval = pow(q0[indexs1] + 
	zA*hparam*f*theta*q1ar[indexs2][ithe-1], f-1.0) 
          - F1mat[ix][ithe];
      yint[ix] = LamF*funcval;

      yint[ix] = F1mat[ix][ithe] + yint[ix]/2.0;

    }

    diffeqsolver(Ncs, yint, qyint);
	
    for(int ix = 0; ix < Nx; ix++) {
      // q0 (s=N)
      indexs1 = ix*(Ns+1) + Ns;
      // q1 (s=2N)
      // q1: theta_index=0:ithe=1
      indexs2 = ix*(2*Ns+1) + 2*Ns;

      funcval = pow(q0[indexs1] + 
	zA*hparam*f*theta*qyint[indexs2], f-1.0) 
          - yint[ix];

      F1new = F1mat[ix][ithe] + LamF*funcval;
  
      // compute error in F1
      errorF += fabs(F1new - F1mat[ix][ithe]);

      F1mat[ix][ithe] = F1new;

      cout<<pow(q0[indexs1] +
              zA*hparam*f*theta*q1ar[indexs2][ithe-1], f-1.0)
              -  F1mat[ix][ithe]<<"  ";
    }

    cout<<endl;

    free(yint);
    free(qyint);

  }
  else if (Fscheme == 5) {

    // Fehlberg method

    double eps = 1e-4;
    double r = 1e6;
    double diff, max;
    double F1new;
    int Ncs = 2*Ns+1;
    double *yint, *qyint, *f1, *f2, *f3;
    yint = (double *) malloc(Nx*sizeof(double));
    qyint = (double *) malloc(Nx*Ncs*sizeof(double));
    f1 = (double *) malloc(Nx*sizeof(double));
    f2 = (double *) malloc(Nx*sizeof(double));
    f3 = (double *) malloc(Nx*sizeof(double));

    while (r > eps) {
     
      for(int ix = 0; ix < Nx; ix++) {
	// q0 (s=N)
	indexs1 = ix*(Ns+1) + Ns;
	// q1 (s=2N)
	// q1: theta_index=0:ithe=1
	indexs2 = ix*(2*Ns+1) + 2*Ns;
	f1[ix] = pow(q0[indexs1] + 
		     zA*hparam*f*theta*q1ar[indexs2][ithe-1], f-1.0) 
          - F1mat[ix][ithe];

	yint[ix] = F1mat[ix][ithe] + hval*f1[ix];

      }

      diffeqsolver(Ncs, yint, qyint);

      for(int ix = 0; ix < Nx; ix++) {
	// q0 (s=N)
	indexs1 = ix*(Ns+1) + Ns;
	// q1 (s=2N)
	// q1: theta_index=0:ithe=1
	indexs2 = ix*(2*Ns+1) + 2*Ns;
	f2[ix] = pow(q0[indexs1] + 
		     zA*hparam*f*theta*qyint[indexs2], f-1.0) 
          - yint[ix];

	yint[ix] = F1mat[ix][ithe] +
	  hval/4.0*(f1[ix]+f2[ix]);

      }

      diffeqsolver(Ncs, yint, qyint);

      for(int ix = 0; ix < Nx; ix++) {
	// q0 (s=N)
	indexs1 = ix*(Ns+1) + Ns;
	// q1 (s=2N)
	// q1: theta_index=0:ithe=1
	indexs2 = ix*(2*Ns+1) + 2*Ns;
	f3[ix] = pow(q0[indexs1] + 
		     zA*hparam*f*theta*qyint[indexs2], f-1.0) 
          - yint[ix];
      }

      max = 1e-18;
      for(int ix = 0; ix < Nx; ix++) {
	diff = fabs( hval/2.0*(f1[ix] + f2[ix]) - 
		     hval/6.0*(f1[ix] + f2[ix] + 4.0*f3[ix]) );
	if (diff > max) max = diff;
      }

      r = max/fabs(hval);
      hval = 0.9 * sqrt(eps/r) * hval;

      cout<<"changing time step to "<<hval<<endl;

    }

    for(int ix = 0; ix < Nx; ix++) {
      F1new = F1mat[ix][ithe] + 
            hval/6.0*(f1[ix] + f2[ix] + 4.0*f3[ix]);
 
      // compute error in F1
      errorF += fabs(F1new - F1mat[ix][ithe]);

      // update F1
      F1mat[ix][ithe] = F1new;

      cout<<pow(q0[indexs1] +
              zA*hparam*f*theta*q1ar[indexs2][ithe-1], f-1.0)
              -  F1mat[ix][ithe]<<"  ";
    }

    cout<<endl;

    free(yint);
    free(qyint);
    free(f1);
    free(f2);
    free(f3);

  }
  else if (Fscheme == 6) {

    // Newton-Raphson

    double tolx, tolf;
    tolx = tolf = MAXERRF;
    int ntrial = MAXITER;  
    int Ncs = 2*Ns+1;
    double temperr;

    double *F1new = (double *) malloc(Nx*sizeof(double));
    double *qF1new = (double *) malloc(Nx*Ncs*sizeof(double));

    // initial guess
    for(int ix = 0; ix < Nx; ix++) {
      F1new[ix] = F1mat[ix][ithe];
    }

    // Newton-Raphson root-finding
    root->mnewt(ntrial, F1new, Nx, tolx, tolf, theta);

    // update q1 with root from Newton-Raphson
    diffeqsolver(Ncs, F1new, qF1new);

    // error
    for(int ix = 0; ix < Nx; ix++) {      
      // q0 (s=N)
      indexs1 = ix*(Ns+1) + Ns;
      // q1 (s=2N)
      // q1: theta_index=0:ithe=1
      indexs2 = ix*(2*Ns+1) + 2*Ns;
      temperr = pow(q0[indexs1] +
              zA*hparam*f*theta*qF1new[indexs2], f-1.0)
	-  F1new[ix];
      cout<<temperr<<"  ";

      // compute error in F1
      errorF += fabs(temperr);

      // update F1
      F1mat[ix][ithe] = F1new[ix];

      // update q1
      for(int ics = 0; ics < Ncs; ics++) {
	q1ar[ix*Ncs + ics][ithe-1] = qF1new[ix*Ncs + ics];
      }

    }

    cout<<endl;
    free(F1new);
    free(qF1new);

  } else if (Fscheme == 7) {

    // Newton-Raphson with line searches and backtracking
    bool check = true;
    int Ncs = 2*Ns+1;
    double temperr;

    double *F1new = (double *) malloc(Nx*sizeof(double));
    double *qF1new = (double *) malloc(Nx*Ncs*sizeof(double));

    // initial guess
    for(int ix = 0; ix < Nx; ix++) {
      F1new[ix] = F1mat[ix][ithe];
    }

    // root-finding
    root->newt(F1new, Nx, check, theta);

    if(check) {
      cout<<"NR with line searches and backtracking failed"<<endl;
      exit(1);
    }

    // update q1 with root
    diffeqsolver(Ncs, F1new, qF1new);

    // error
    for(int ix = 0; ix < Nx; ix++) {      
      // q0 (s=N)
      indexs1 = ix*(Ns+1) + Ns;
      // q1 (s=2N)
      // q1: theta_index=0:ithe=1
      indexs2 = ix*(2*Ns+1) + 2*Ns;
      temperr = pow(q0[indexs1] +
              zA*hparam*f*theta*qF1new[indexs2], f-1.0)
	-  F1new[ix];
      cout<<temperr<<"  ";

      // compute error in F1
      errorF += fabs(temperr);

      // update F1
      F1mat[ix][ithe] = F1new[ix];

      // update q1
      for(int ics = 0; ics < Ncs; ics++) {
	q1ar[ix*Ncs + ics][ithe-1] = qF1new[ix*Ncs + ics];
      }

    }

    cout<<endl;
    free(F1new);
    free(qF1new);

  } else if (Fscheme == 8) {

    // Broyden's method
    bool check = true;
    int Ncs = 2*Ns+1;
    double temperr;

    double *F1new = (double *) malloc(Nx*sizeof(double));
    double *qF1new = (double *) malloc(Nx*Ncs*sizeof(double));

    // initial guess
    for(int ix = 0; ix < Nx; ix++) {
      F1new[ix] = F1mat[ix][ithe];
    }

    // root-finding
    root->broydn(F1new, Nx, check, theta);

    if(check) {
      cout<<"Broyden's method failed"<<endl;
      exit(1);
    }

    // update q1 with root
    diffeqsolver(Ncs, F1new, qF1new);

    // error
    for(int ix = 0; ix < Nx; ix++) {      
      // q0 (s=N)
      indexs1 = ix*(Ns+1) + Ns;
      // q1 (s=2N)
      // q1: theta_index=0:ithe=1
      indexs2 = ix*(2*Ns+1) + 2*Ns;
      temperr = pow(q0[indexs1] +
              zA*hparam*f*theta*qF1new[indexs2], f-1.0)
	-  F1new[ix];
      cout<<temperr<<"  ";

      // compute error in F1
      errorF += fabs(temperr);

      // update F1
      F1mat[ix][ithe] = F1new[ix];

      // update q1
      for(int ics = 0; ics < Ncs; ics++) {
	q1ar[ix*Ncs + ics][ithe-1] = qF1new[ix*Ncs + ics];
      }

    }

    cout<<endl;
    free(F1new);
    free(qF1new);

  }   
  else {
    double F1new;
    for(int ix = 0; ix < Nx; ix++) {
      // q0 (s=N)
      indexs1 = ix*(Ns+1) + Ns;
      // q1 (s=2N)
      // q1: theta_index=0:ithe=1
      indexs2 = ix*(2*Ns+1) + 2*Ns;
      if(Fscheme == 0) {
        // explicit
        F1new = (1.0 - LamF)*F1mat[ix][ithe] +
	  LamF * pow(q0[indexs1] + 
           zA*hparam*f*theta*q1ar[indexs2][ithe-1], f-1.0); 
      } else if (Fscheme == 1) {
        // semi-implicit
        F1new = F1mat[ix][ithe] + 
	   LamF*pow(q0[indexs1] + 
              zA*hparam*f*theta*q1ar[indexs2][ithe-1], f-1.0);
        F1new /= (1.0 + LamF);
      } else {
	cout<<"unknown Fscheme"<<endl;
	exit(1);
      }

      // compute error in F1
      errorF += fabs(F1new - F1mat[ix][ithe]);

      cout<<pow(q0[indexs1] +
              zA*hparam*f*theta*q1ar[indexs2][ithe-1], f-1.0)
              -  F1mat[ix][ithe]<<"  ";

       
      // update F1
      F1mat[ix][ithe] = F1new;
    }    

     cout<<endl;
  }

  // average error
  errorF /= double(Nx);
  return errorF;
}


void Confinpg::calcphi() 
{
  int index_11, index_10;      
  double tF1f, tF0f;
  double integ_s;
  double tF1f0, tF0fN;
  // theta=1
  int itheta = Ntheta - 1;

  // integrate over s for each x

  for(int ix = 0; ix < Nx; ix++) {

    // integrate over s
    integ_s = 0.0;
    for(int ics = 0; ics <= Ns; ics++) {
      tF1f = q0[ix*(Ns+1)+ics];
      // 1+s
      index_11 = ix*(2*Ns+1)+Ns+ics;
      tF1f += zA*hparam*f*q1t[index_11][itheta-1];
   
      // 1-s
      index_10 = ix*(2*Ns+1)+Ns-ics;
      tF0f = q1t[index_10][itheta-1];

      if (ics == 0) {
	tF0fN = tF0f;
        tF1f0 = tF1f;
      }

      if (ics==0 || ics==Ns) 
        integ_s +=  3.0/8.0*tF1f*tF0f;
      else if (ics==1 || ics==Ns-1)
	integ_s +=  7.0/6.0*tF1f*tF0f;
      else if (ics==2 || ics==Ns-2)
	integ_s +=  23.0/24.0*tF1f*tF0f;
      else  
	integ_s +=  tF1f*tF0f;

      // end for s-loop
    }

    integ_s *= Dels;

    // phitot, total end densities
    phitot[ix] = zA*f*integ_s;
    phie[ix] = zA*f/NA*tF1f0*tF0fN;
    phiue[ix] = zA*f/NA*tF0fN;

    // end for ix-loop
  }
}


void Confinpg::calcphisol() 
{
  int index_11, index_10;      
  double tF1f, tF0f;
  double integ_s;
  double tF1f0, tF0fN;
  // theta=1
  int itheta = Ntheta - 1;

  // integrate over s for each x

  for(int ix = 0; ix < Nx; ix++) {

    // integrate over s
    integ_s = 0.0;
    for(int ics = 0; ics <= Ns; ics++) {
      tF1f = q0[ix*(Ns+1)+ics];
      // 1+s
      index_11 = ix*(2*Ns+1)+Ns+ics;
      tF1f += zA*hparam*f*q1s[index_11][itheta-1];
   
      // 1-s
      index_10 = ix*(2*Ns+1)+Ns-ics;
      tF0f = q1s[index_10][itheta-1];

      if (ics == 0) {
	tF0fN = tF0f;
        tF1f0 = tF1f;
      }

      if (ics==0 || ics==Ns) 
        integ_s +=  3.0/8.0*tF1f*tF0f;
      else if (ics==1 || ics==Ns-1)
	integ_s +=  7.0/6.0*tF1f*tF0f;
      else if (ics==2 || ics==Ns-2)
	integ_s +=  23.0/24.0*tF1f*tF0f;
      else  
	integ_s +=  tF1f*tF0f;

      // end for s-loop
    }

    integ_s *= Dels;

    // phisol, sol end densities
    phisol[ix] = zA*f*integ_s;
    phiesol[ix] = zA*f/NA*tF1f0*tF0fN;
    phiuesol[ix] = zA*f/NA*tF0fN;

    // end for ix-loop
  }
}


double Confinpg::calcavphi()
{

  double avphi = 0.0;

  // calculate fourier sine coeffs of phitot

  double *xdata;
  xdata = (double *) fftw_malloc(sizeof(double)*Nx);

  for(int ix = 0; ix < Nx; ix++) {
    xdata[ix] = phitot[ix];
  }

  fftw_execute_r2r(fftwp, xdata, xdata);


  // xdata contains un-normalized fourier coeffs
  // to normalize, divide by 2(Nx+1)

  // calculate average density
  // loop over even values of index m
  for(int m = 0; m < Nx; m += 2) {
    avphi += xdata[m]/(double(m) + 1.0);
  }
  avphi *= 2.0/PI/(double(Nx)+1.0);

  fftw_free(xdata);

  return avphi;
}


void Confinpg::diffsolve_0()
{
  // initial condition=1.0, s in [0,1]
  int Ncs = Ns + 1;
  double *qres = (double *) malloc(Nx*Ncs*sizeof(double));
  double *qinit;
  qinit = (double *) malloc(Nx*sizeof(double));
  for(int ix = 0; ix < Nx; ix++) {
    qinit[ix] = 1.0;
  }

  diffeqsolver(Ncs, qinit, qres);

  for(int i = 0; i < Nx*Ncs; i++) {
    q0[i] = qres[i];
  }

  free(qres);
  free(qinit);
}


void Confinpg::diffsolve_1s(int itheta)
{
  // initial condition=F1s, s in [0,2]
  int Ncs = 2*Ns + 1;
  double *qres = (double *) malloc(Nx*Ncs*sizeof(double));
  double *qinit;
  qinit = (double *) malloc(Nx*sizeof(double));
  for(int ix = 0; ix < Nx; ix++) {
    qinit[ix] = F1s[ix][itheta];
  }

  diffeqsolver(Ncs, qinit, qres);

  for(int i = 0; i < Nx*Ncs; i++) {
    q1s[i][itheta-1] = qres[i];
  }

  free(qres);
  free(qinit);
}


void Confinpg::diffsolve_1t(int itheta)
{
  // initial condition=F1t, s in [0,2]
  int Ncs = 2*Ns + 1;
  double *qres = (double *) malloc(Nx*Ncs*sizeof(double));
  double *qinit;
  qinit = (double *) malloc(Nx*sizeof(double));
  for(int ix = 0; ix < Nx; ix++) {
    qinit[ix] = F1t[ix][itheta];
  }

  diffeqsolver(Ncs, qinit, qres);

  for(int i = 0; i < Nx*Ncs; i++) {
    q1t[i][itheta-1] = qres[i];
  }

  free(qres);
  free(qinit);
}


void Confinpg::diffsolve_1t0()
{
  // initial condition=q_0(s=1)^{f-1}=F1s,t at theta=0, s in [0,1]
  int Ncs = Ns + 1;
  double *qres = (double *) malloc(Nx*Ncs*sizeof(double));
  double *qinit;
  qinit = (double *) malloc(Nx*sizeof(double));
  for(int ix = 0; ix < Nx; ix++) {
    qinit[ix] = F1s[ix][0];
  }

  diffeqsolver(Ncs, qinit, qres);

  for(int i = 0; i < Nx*Ncs; i++) {
    q1t0[i] = qres[i];
  }

  free(qres);
  free(qinit);
}


void Confinpg::diffeqsolver(int Nssize, double *qic, double *qres)
{
  // fourier coeffs (n+1/3)
  double *qint = (double *) fftw_malloc(Nx*sizeof(double));
  // initial condition, s=0
  for(int ix = 0; ix < Nx; ix++){
    qres[ix*Nssize] = qic[ix];
  }

  // step along contour
  for (int ics = 1; ics < Nssize; ics++) {

    // 1/3rd step
    // real space
    // using previous contour position as starting point
    for(int ix = 0; ix < Nx; ix++) {
      qint[ix] = exp(-Dels/2.0 * wA[ix]) * qres[ix*Nssize+ics-1];
    }

    // fourier coeffs (n+1/3)
    // DST-I
    fftw_execute_r2r(fftwp, qint, qint);

    // normalization
    // fourier coeffs (n+2/3)
    for(int m = 0; m < Nx; m++) {
      qint[m] /= (2.0*(double(Nx) + 1.0));
      qint[m] *= exp(-Dels * PI*PI*double(m+1)*double(m+1)/len/len );
    }

    // real (n+2/3)
    // DST-I
    fftw_execute_r2r(fftwp, qint, qint);

    // real (n+1)
    for(int ix = 0; ix < Nx; ix++) {
      qres[ix*Nssize + ics] = exp(-Dels/2.0 * wA[ix]) * qint[ix];
    }
  }
  
  fftw_free(qint);
}


void Confinpg::calcalpha() 
{
  for(int ix = 0; ix < Nx; ix++) {
    alpha[ix] = 1.0 - phiue[ix]/phie[ix];
  }
}


void Confinpg::calcalphasol() 
{
  for(int ix = 0; ix < Nx; ix++) {
    alphasol[ix] = 1.0 - phiuesol[ix]/phiesol[ix];
  }
}


void Confinpg::calcmoments()
{
  calcF0t();

  // total moments 
  // number-average mol wt

  double term1, term2;
  double *F0xvals = (double *) fftw_malloc(Nx*sizeof(double));

  // theta=1  
  for(int ix = 0; ix < Nx; ix++) {
    F0xvals[ix] = F0t[ix][Ntheta-1];
  }

  // fourier coeffs  
  fftw_execute_r2r(fftwp, F0xvals, F0xvals);

  // integral over x, theta=1
  term1 = 0.0;
  for(int m = 0; m < Nx; m += 2) {
    term1 += F0xvals[m]/(double(m) + 1.0);
  }

  // integral over theta
  for(int ix = 0; ix < Nx; ix++) {
    F0xvals[ix] = 0.0;
    for(int itheta = 0; itheta < Ntheta; itheta++) {
      if (itheta==0 || itheta==Ntheta-1) 
	F0xvals[ix] +=  3.0/8.0*F0t[ix][itheta];
      else if (itheta==1 || itheta==Ntheta-2)
	F0xvals[ix] +=  7.0/6.0*F0t[ix][itheta];
      else if (itheta==2 || itheta==Ntheta-3)
	F0xvals[ix] +=  23.0/24.0*F0t[ix][itheta];
      else  
	F0xvals[ix] +=  F0t[ix][itheta];
    }
    F0xvals[ix] *= DelTheta;
  }

  // integral over x
  // fourier coeffs  
  fftw_execute_r2r(fftwp, F0xvals, F0xvals);  
  term2 = 0.0;
  for(int m = 0; m < Nx; m += 2) {
    term2 += F0xvals[m]/(double(m) + 1.0);
  }    

  Nn = term1/term2;

  // weight-average mol wt

  double term3;
  // derivative wrt theta at theta=1
  for(int ix = 0; ix < Nx; ix++) {
    F0xvals[ix] = (3.0*F0t[ix][Ntheta-1] - 4.0*F0t[ix][Ntheta-2] 
            + F0t[ix][Ntheta-3])/2.0/DelTheta;
  }

  // integral over x
  // fourier coeffs  
  fftw_execute_r2r(fftwp, F0xvals, F0xvals);  
  term3 = 0.0;
  for(int m = 0; m < Nx; m += 2) {
    term3 += F0xvals[m]/(double(m) + 1.0);
  }    
  
  Nw = term3/term1 + 1.0;

  fftw_free(F0xvals);

}


void Confinpg::calcmomentsol()
{
  calcF0s();

  // moments in sol phase
  // number-average mol wt

  double term1, term2;
  double *F0xvals = (double *) fftw_malloc(Nx*sizeof(double));

  // theta=1  
  for(int ix = 0; ix < Nx; ix++) {
    F0xvals[ix] = F0s[ix][Ntheta-1];
  }
  // fourier coeffs  
  fftw_execute_r2r(fftwp, F0xvals, F0xvals);

  // integral over x, theta=1
  term1 = 0.0;
  for(int m = 0; m < Nx; m += 2) {
    term1 += F0xvals[m]/(double(m) + 1.0);
  }

  // integral over theta
  for(int ix = 0; ix < Nx; ix++) {
    F0xvals[ix] = 0.0;
    for(int itheta = 0; itheta < Ntheta; itheta++) {
      if (itheta==0 || itheta==Ntheta-1) 
	F0xvals[ix] +=  3.0/8.0*F0s[ix][itheta];
      else if (itheta==1 || itheta==Ntheta-2)
	F0xvals[ix] +=  7.0/6.0*F0s[ix][itheta];
      else if (itheta==2 || itheta==Ntheta-3)
	F0xvals[ix] +=  23.0/24.0*F0s[ix][itheta];
      else  
	F0xvals[ix] +=  F0s[ix][itheta];
    }
    F0xvals[ix] *= DelTheta;
  }
  // integral over x
  // fourier coeffs  
  fftw_execute_r2r(fftwp, F0xvals, F0xvals);  
  term2 = 0.0;
  for(int m = 0; m < Nx; m += 2) {
    term2 += F0xvals[m]/(double(m) + 1.0);
  }    

  Nnsol = term1/term2;

  // weight-average mol wt

  double term3;
  // derivative wrt theta at theta=1
  for(int ix = 0; ix < Nx; ix++) {
    F0xvals[ix] = (3.0*F0s[ix][Ntheta-1] - 4.0*F0s[ix][Ntheta-2] 
            + F0s[ix][Ntheta-3])/2.0/DelTheta;
  }
  // integral over x
  // fourier coeffs  
  fftw_execute_r2r(fftwp, F0xvals, F0xvals);  
  term3 = 0.0;
  for(int m = 0; m < Nx; m += 2) {
    term3 += F0xvals[m]/(double(m) + 1.0);
  }    
  
  Nwsol = term3/term1 + 1.0;

  fftw_free(F0xvals);

}


void Confinpg::calcF0s()
{
  for(int ix = 0; ix < Nx; ix++) {
    for(int itheta = 0; itheta < Ntheta; itheta++) {
      F0s[ix][itheta] = pow(F1s[ix][itheta], f/(f-1.0));
    }
  }
}


void Confinpg::calcF0t()
{
  for(int ix = 0; ix < Nx; ix++) {
    for(int itheta = 0; itheta < Ntheta; itheta++) {
      F0t[ix][itheta] = pow(F1t[ix][itheta], f/(f-1.0));
    }
  }
}


void Confinpg::calcphiumon()
{
  double integral, integrand;
  // calculate q1t0 
  diffsolve_1t0();

  // integrate over s for each x

  for(int ix = 0; ix < Nx; ix++) {

    // integrate over s
    integral = 0.0;
    for(int ics = 0; ics <= Ns; ics++) {
      integrand = q0[ix*(Ns+1)+Ns-ics]*q1t0[ix*(Ns+1)+ics];
      if (ics==0 || ics==Ns) 
        integral += 3.0/8.0*integrand;
      else if (ics==1 || ics==Ns-1)
	integral += 7.0/6.0*integrand;
      else if (ics==2 || ics==Ns-2)
	integral += 23.0/24.0*integrand;
      else  
	integral += integrand;
 
      // end for s-loop
    }

    integral *= Dels;

    // unreacted monomer density
    phiunrmon[ix] = zA*f*integral;

    // end for ix-loop
  }

}


void Confinpg::printres() 
{
  calcalpha();
  calcphisol();
  calcalphasol();
  calcphiumon();
  calcmoments();
  calcmomentsol();

  double xval;
  ofstream fout;
  fout.open("outfile.dat", ios::binary|ios::out);

  for(int ix = 0; ix < Nx; ix++) {
    xval = len*(double(ix)+1.0)/(double(Nx)+1.0);
    fout<<xval<<"  "<<phitot[ix]<<"  "<<phisol[ix]<<"  "
	<<alpha[ix]<<"  "<<alphasol[ix]<<"  "
        <<phiunrmon[ix]<<"  "<<wA[ix]<<endl; 
  }
   
  fout.close();

  fout.open("wAin.dat", ios::trunc|ios::out);
  for(int ix = 0; ix < Nx; ix++) {
    fout<<wA[ix]<<endl; 
  }
  fout.close();

  // F1s,t at theta=1
  fout.open("F1vals.dat", ios::binary|ios::out);
  for(int ix = 0; ix < Nx; ix++) {
    xval = len*(double(ix)+1.0)/(double(Nx)+1.0);
    fout<<xval<<"  "<<F1s[ix][Ntheta-1]<<"  "
	<<F1t[ix][Ntheta-1]<<endl; 
  }
   
  fout.close();

  cout<<"Nn = "<<Nn<<endl;
  cout<<"Nw = "<<Nw<<endl;
  cout<<"Nnsol = "<<Nnsol<<endl;
  cout<<"Nwsol = "<<Nwsol<<endl;
  cout<<"zA = "<<zA<<endl;

}


void Confinpg::calcfns(double *F1, double *funcF1, double theta) 
{
  int Ncs = 2*Ns+1;
  double *qF1 = (double *) malloc(Nx*Ncs*sizeof(double));
  int indexs1, indexs2;

  diffeqsolver(Ncs, F1, qF1);

  for(int ix = 0; ix < Nx; ix++) {
      // q0 (s=N)
      indexs1 = ix*(Ns+1) + Ns;
      // q1 (s=2N)
      // q1: theta_index=0:ithe=1
      indexs2 = ix*(2*Ns+1) + 2*Ns;
      funcF1[ix] = pow(q0[indexs1] + 
              zA*hparam*f*theta*qF1[indexs2], f-1.0) - F1[ix];
  }

  free(qF1);
}


void Confinpg::calcallfns(double *xvec, double *fvec) 
{

  // 0 to Nx-1: F1t
  // Nx to 2*Nx-1: w
  // 2*Nx: z

  double errw = 0.0;
  int indexs1, indexs2;

  for(int ix = 0; ix < Nx; ix++) {
    F1t[ix][Ntheta-1] = xvec[ix];
    wA[ix] = xvec[Nx+ix];
  }
  zA = xvec[2*Nx];

  diffsolve_0();
  diffsolve_1t(Ntheta-1);

  for(int ix = 0; ix < Nx; ix++) {
      // q0 (s=N)
      indexs1 = ix*(Ns+1) + Ns;
      // q1 (s=2N)
      indexs2 = ix*(2*Ns+1) + 2*Ns;
      fvec[ix] = pow(q0[indexs1] + 
              zA*hparam*f*q1t[indexs2][Ntheta-2], f-1.0) 
               - F1t[ix][Ntheta-1];
  }

  calcphi();

  for(int ix = 0; ix < Nx; ix++) {
    fvec[Nx+ix] = -wA[ix]/uparam + phitot[ix];
    errw +=  fvec[Nx+ix];
  }
  errw /= double(Nx);

  double avphi = calcavphi();
  fvec[2*Nx] = 1.0 - avphi;

  cout<<errw<<"  "<<1.0-avphi<<endl;

}
