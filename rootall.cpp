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
#include "rootall.h"

using namespace std;

RootAll::RootAll(ConfGel *cgel) :Pointers(cgel) {}

void RootAll::lubksb(double **a, int n, int *indx, double *b)
{
  int i,ii=0,ip,j;
  double sum;

  for (i=0;i<n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii!=0)
      for (j=ii-1;j<i;j++) sum -= a[i][j]*b[j];
    else if (sum!=0) ii=i+1;
    b[i]=sum;
  }
  for (i=n-1;i>=0;i--) {
    sum=b[i];
    for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}


void RootAll::ludcmp(double **a, int n, int *indx, double &d)
{
  const double TINY = 1.0e-20;
  int i,imax,j,k;
  double big,dum,sum,temp;
  double *vv;

  vv = (double *) malloc(n*sizeof(double));

  d=1.0;
  for (i=0;i<n;i++) {
    big=0.0;
    for (j=0;j<n;j++)
      if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big == 0.0) {
      cout<<"Singular matrix in routine ludcmp"<<endl;
      free(vv);
      exit(1);
    }

    vv[i]=1.0/big;
  }
  for (j=0;j<n;j++) {
    for (i=0;i<j;i++) {
      sum=a[i][j];
      for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<n;i++) {
      sum=a[i][j];
      for (k=0;k<j;k++)
	sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {
      for (k=0;k<n;k++) {
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      d = -d;
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    if (j != n-1) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<n;i++) a[i][j] *= dum;
    }
  }
  free(vv);
}


void RootAll::mnewt(int ntrial, double *x, int n, 
		    double tolx, double tolf)
{
  int k,i;
  double errx,errf,d,*fvec,**fjac,*p;

  int *indx=(int *) malloc(n*sizeof(int));
  p = (double *) malloc(n*sizeof(double));
  fvec=(double *) malloc(n*sizeof(double));
  fjac=mem->create_2d_double_array(n,n);

  for (k=0;k<ntrial;k++) {
    usrfun(x,n,fvec,fjac);
    errf=0.0;
    for (i=0;i<n;i++) errf += fabs(fvec[i]);
    if (errf <= tolf) {
      free(indx);
      free(p);
      free(fvec);
      mem->destroy_2d_double_array(fjac);
      cout<<k<<" Newton iterations"<<endl;
      return;		 
    }
                      
    for (i=0;i<n;i++) p[i] = -fvec[i];
    ludcmp(fjac,n,indx,d);
    lubksb(fjac,n,indx,p);
    errx=0.0;
    for (i=0;i<n;i++) {
      errx += fabs(p[i]);
      x[i] += p[i];
    }
    if (errx <= tolx) {
      free(indx);
      free(p);
      free(fvec);
      mem->destroy_2d_double_array(fjac);
      cout<<k<<" Newton iterations"<<endl;
      return;
    }
  }
  free(indx);
  free(p);
  free(fvec);
  mem->destroy_2d_double_array(fjac);
  cout<<"Newton-Raphson failed to converge"<<endl;
  exit(1);
}


void RootAll::fdjac(int n, double *x, double *fvec, double **df)
{

  const double EPS = 1.0e-4;
  int i,j;
  double h,temp,*f;
  f = (double *) malloc(n*sizeof(double));

  for (j=0;j<n;j++) {
    temp=x[j];
    h=EPS*fabs(temp);
    if (h == 0.0) h=EPS;
    x[j]=temp+h;
    h=x[j]-temp;
    confpg->calcallfns(x,f);
    x[j]=temp;
    for (i=0;i<n;i++) df[i][j]=(f[i]-fvec[i])/h;
  }

  free(f);
}


void RootAll::usrfun(double *x, int n, double *fx, double **fjac)
{
  
  confpg->calcallfns(x, fx);

  fdjac(n, x, fx, fjac);

}


void RootAll::lnsrch(int n, double *xold, double fold, 
                     double *g, double *p, double *x, 
                     double &f, const double stpmax, bool &check,
                     double *fvec)
{

  const double ALF = 1.0e-4;
  const double TOLX = 1.0e-8;
  int i;
  double a,alam,alam2=0.0,alamin,b,disc,f2=0.0;
  double rhs1,rhs2,slope,sum,temp,test,tmplam;

  check=false;
  sum=0.0;
  for (i=0;i<n;i++) sum += p[i]*p[i];
  sum=sqrt(sum);
  if (sum > stpmax)
    for (i=0;i<n;i++) p[i] *= stpmax/sum;
  slope=0.0;
  for (i=0;i<n;i++)
    slope += g[i]*p[i];
  if(slope >= 0.0) {
    cout<<"roundoff problem in lnsrch"<<endl;
    exit(1);
  }
  test=0.0;
  for (i=0;i<n;i++) {
    temp=fabs(p[i])/MAX(fabs(xold[i]),1.0);
    if (temp > test) test=temp;
  }
  alamin=TOLX/test;
  alam=1.0;
  for (;;) {
    for (i=0;i<n;i++) x[i]=xold[i]+alam*p[i];
    f=fmin(x,fvec,n);
    if (alam < alamin) {
      for (i=0;i<n;i++) x[i]=xold[i];
      check=true;
      return;
    } else if (f <= fold+ALF*alam*slope) return;
    else {
      if (alam == 1.0)
	tmplam = -slope/(2.0*(f-fold-slope));
      else {
	rhs1 = f-fold-alam*slope;
	rhs2=f2-fold-alam2*slope;
	a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/
	  (alam-alam2);
	if (a == 0.0) tmplam = -slope/(2.0*b);
	else {
	  disc=b*b-3.0*a*slope;
	  if (disc<0.0) tmplam = 0.5*alam;
	  else if (b<=0.0) tmplam = (-b+sqrt(disc))/(3.0*a);
	  else tmplam=-slope/(b+sqrt(disc));
	}
	if (tmplam>0.5*alam)
	  tmplam=0.5*alam;
      }
    }
    alam2=alam;
    f2 = f;
    alam=MAX(tmplam,0.1*alam);
  }
}


double RootAll::MAX(double a, double b) 
{
  return ( (a>b)? a:b);
}


void RootAll::newt(double *x, int n, bool &check)
{
  const int MAXITS = 200;
  const double TOLF=1.0e-8, TOLMIN=1.0e-12, STPMX=100.0;
  const double TOLX = 1.0e-8;

  int i,its,j; 
  double d,den,f,fold,stpmax,sum,temp,test;

  int *indx = (int *) malloc(n*sizeof(int));
  double *g = (double *) malloc(n*sizeof(double));
  double *p = (double *) malloc(n*sizeof(double));
  double *xold = (double *) malloc(n*sizeof(double));

  double **fjac = mem->create_2d_double_array(n,n);
  double *fvec =  (double *) malloc(n*sizeof(double));
  f = fmin(x,fvec,n);
  test=0.0;
  for (i=0;i<n;i++)
    if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
  if (test<0.01*TOLF){
    check=false;
    free(indx);
    free(g);
    free(p);
    free(xold);
    mem->destroy_2d_double_array(fjac);
    free(fvec);
    return;
  }
  sum=0.0;
  for (i=0;i<n;i++) sum += x[i]*x[i];
  stpmax=STPMX*MAX(sqrt(sum),(double)n);
  for (its=0;its<MAXITS;its++) {

    cout<<"iteration # "<<its<<endl;

    fdjac(n,x,fvec,fjac);
    for (i=0;i<n;i++) {
      sum=0.0;
      for (j=0;j<n;j++) sum += fjac[j][i]*fvec[j];
      g[i]=sum;
    }
    for (i=0;i<n;i++) xold[i]=x[i];
    fold=f;
    for (i=0;i<n;i++) p[i] = -fvec[i];
    ludcmp(fjac,n,indx,d);
    lubksb(fjac,n,indx,p);
    lnsrch(n,xold,fold,g,p,x,f,stpmax,check,fvec);
    test=0.0;
    for (i=0;i<n;i++)
      if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
    if (test < TOLF) {
      check=false;    
      free(indx);
      free(g);
      free(p);
      free(xold);
      mem->destroy_2d_double_array(fjac);
      free(fvec);
      return;
    }
    if (check) {
      test=0.0;
      den=MAX(f,0.5*n);
      for (i=0;i<n;i++) {
	temp=fabs(g[i])*MAX(fabs(x[i]),1.0)/den;
	if (temp > test) test=temp;
      }
      check=(test < TOLMIN);
      free(indx);
      free(g);
      free(p);
      free(xold);
      mem->destroy_2d_double_array(fjac);
      free(fvec);
      return;
    }
    test=0.0;
    for (i=0;i<n;i++) {
      temp=(fabs(x[i]-xold[i]))/MAX(fabs(x[i]),1.0);
      if (temp > test) test=temp;
    }
    if (test < TOLX) {    
      free(indx);
      free(g);
      free(p);
      free(xold);
      mem->destroy_2d_double_array(fjac);
      free(fvec);
      return;
    }
  }
  cout<<"MAXITS exceeded in newt"<<endl;
  free(indx);
  free(g);
  free(p);
  free(xold);
  mem->destroy_2d_double_array(fjac);
  free(fvec);
  exit(1);
}


double RootAll::fmin(double *x, double *fx, int n)
{
  double sum;
  confpg->calcallfns(x, fx);
  sum = 0.0;
  for (int i=0;i<n;i++) sum += fx[i]*fx[i];
  return 0.5*sum;
}


void RootAll::qrupdt(double **r, double **qt, int n, double *u, double *v)
{
  int i,k;

  for (k=n-1;k>=0;k--) {
    if (u[k]!=0.0) break;
  }
  if (k < 0) k=0;
  for (i=k-1;i>=0;i--) {
    rotate(r,qt,n,i,u[i],-u[i+1]);
    if (u[i] == 0.0) u[i]=fabs(u[i+1]);
    else if (fabs(u[i]) > fabs(u[i+1]))
      u[i]=fabs(u[i])*sqrt(1.0+(u[i+1]/u[i])*(u[i+1]/u[i]));
    else u[i]=fabs(u[i+1])*sqrt(1.0+(u[i]/u[i+1])*(u[i]/u[i+1]));
  }
  for (i=0;i<n;i++) r[0][i] += u[0]*v[i];
  for (i=0;i<k;i++)
    rotate(r,qt,n,i,r[i][i],-r[i+1][i]);
}


void RootAll::qrdcmp(double **a, int n, double *c, double *d, bool &sing)
{
  int i,j,k;
  double scale,sigma,sum,tau;

  sing=false;
  for (k=0;k<n-1;k++) {
    scale=0.0;
    for (i=k;i<n;i++) scale=MAX(scale,fabs(a[i][k]));
    if (scale == 0.0) {
      sing=true;
      c[k]=d[k]=0.0;
    } else {
      for (i=k;i<n;i++) a[i][k] /= scale;
      for (sum=0.0,i=k;i<n;i++) sum += a[i][k]*a[i][k];
      sigma=SIGN(sqrt(sum),a[k][k]);
      a[k][k] += sigma;
      c[k]=sigma*a[k][k];
      d[k] = -scale*sigma;
      for (j=k+1;j<n;j++) {
	for (sum=0.0,i=k;i<n;i++) sum += a[i][k]*a[i][j];
	tau=sum/c[k];
	for (i=k;i<n;i++) a[i][j] -= tau*a[i][k];
      }
    }
  }
  d[n-1]=a[n-1][n-1];
  if (d[n-1] == 0.0) sing=true;
}


void RootAll::rotate(double **r, double **qt, int n, int i, double a, double b)
{
  int j;
  double c,fact,s,w,y;

  if (a == 0.0) {
    c=0.0;
    s=(b >= 0.0 ? 1.0 : -1.0);
  } else if (fabs(a) > fabs(b)) {
    fact=b/a;
    c=SIGN(1.0/sqrt(1.0+(fact*fact)),a);
    s=fact*c;
  } else {
    fact=a/b;
    s=SIGN(1.0/sqrt(1.0+(fact*fact)),b);
    c=fact*s;
  }
  for (j=i;j<n;j++) {
    y=r[i][j];
    w=r[i+1][j];
    r[i][j]=c*y-s*w;
    r[i+1][j]=s*y+c*w;
  }
  for (j=0;j<n;j++) {
    y=qt[i][j];
    w=qt[i+1][j];
    qt[i][j]=c*y-s*w;
    qt[i+1][j]=s*y+c*w;
  }
}


void RootAll::rsolv(double **a, int n, double *d, double *b)
{
  int i,j;
  double sum;

  b[n-1] /= d[n-1];
  for (i=n-2;i>=0;i--) {
    for (sum=0.0,j=i+1;j<n;j++) sum += a[i][j]*b[j];
    b[i]=(b[i]-sum)/d[i];
  }
}


double RootAll::SIGN(double a, double b) 
{
  return ((b) >= 0.0 ? fabs(a) : -fabs(a));
}


void RootAll::broydn(double *x, int n, bool &check) 
{
  const int MAXITS=200;
  const double EPS=1.0e-8;
  const double TOLF=1.0e-8;
  const double TOLX=EPS;
  const double STPMX=100.0;
  const double TOLMIN=1.0e-12;
  bool restrt, sing, skip;
  int i,its,j,k;
  double den,f,fold,stpmax,sum,temp,test;

  double **qt, **r;
  qt = mem->create_2d_double_array(n,n);
  r = mem->create_2d_double_array(n,n);
  double *c, *d, *fvcold, *g, *p, *s, *t, *w, *xold;
  c = (double *) malloc(n*sizeof(double));
  d = (double *) malloc(n*sizeof(double));
  fvcold = (double *) malloc(n*sizeof(double));
  g = (double *) malloc(n*sizeof(double));
  p = (double *) malloc(n*sizeof(double));
  s = (double *) malloc(n*sizeof(double));
  t = (double *) malloc(n*sizeof(double));
  w = (double *) malloc(n*sizeof(double));
  xold = (double *) malloc(n*sizeof(double));

  double *fvec = (double *) malloc(n*sizeof(double));

  f=fmin(x,fvec,n);
  test=0.0;
  for (i=0;i<n;i++)
    if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
  if (test<0.01*TOLF) {
    check=false;
    mem->destroy_2d_double_array(qt);
    mem->destroy_2d_double_array(r);
    free(c);
    free(d);
    free(fvcold);
    free(g);
    free(p);
    free(s);
    free(t);
    free(w);
    free(xold);
    free(fvec);
    return;
  }

  for (sum=0.0,i=0;i<n;i++) sum += x[i]*x[i];
  stpmax=STPMX*MAX(sqrt(sum),(double) n );
  restrt=true;
  for (its=1;its<=MAXITS;its++) {

    cout<<"iteration # "<<its<<endl;

    if (restrt) {
      fdjac(n,x,fvec,r);
      qrdcmp(r,n,c,d,sing);
      if (sing) {
	cout<<"singular Jacobian in broydn"<<endl;
	mem->destroy_2d_double_array(qt);
	mem->destroy_2d_double_array(r);
	free(c);
	free(d);
	free(fvcold);
	free(g);
	free(p);
	free(s);
	free(t);
	free(w);
	free(xold);
	free(fvec);
	return;
      }

      for (i=0;i<n;i++) {
	for (j=0;j<n;j++) qt[i][j]=0.0;
	qt[i][i]=1.0;
      }
      for (k=0;k<n-1;k++) {
	if (c[k] != 0.0) {
	  for (j=0;j<n;j++) {
	    sum=0.0;
	    for (i=k;i<n;i++)
	      sum += r[i][k]*qt[i][j];
	    sum /= c[k];
	    for (i=k;i<n;i++)
	      qt[i][j] -= sum*r[i][k];
	  }
	}
      }
      for (i=0;i<n;i++) {
	r[i][i]=d[i];
	for (j=0;j<i;j++) r[i][j]=0.0;
      }
    } else {
      for (i=0;i<n;i++) s[i]=x[i]-xold[i];
      for (i=0;i<n;i++) {
	for (sum=0.0,j=i;j<n;j++) sum += r[i][j]*s[j];
	t[i]=sum;
      }
      skip=true;
      for (i=0;i<n;i++) {
	for (sum=0.0,j=0;j<n;j++) sum += qt[j][i]*t[j];
	w[i]=fvec[i]-fvcold[i]-sum;
	if (fabs(w[i]) >= EPS*(fabs(fvec[i])+fabs(fvcold[i]))) 
	  skip=false;
	else w[i]=0.0;
      }
      if (!skip) {
	for (i=0;i<n;i++) {
	  for (sum=0.0,j=0;j<n;j++) sum += qt[i][j]*w[j];
	  t[i]=sum;
	}
	for (den=0.0,i=0;i<n;i++) den += s[i]*s[i];
	for (i=0;i<n;i++) s[i] /= den;
	qrupdt(r,qt,n,t,s);
	for (i=0;i<n;i++) {
	  if (r[i][i] == 0.0) {
	    cout<<"r singular in broydn"<<endl;
	    mem->destroy_2d_double_array(qt);
	    mem->destroy_2d_double_array(r);
	    free(c);
	    free(d);
	    free(fvcold);
	    free(g);
	    free(p);
	    free(s);
	    free(t);
	    free(w);
	    free(xold);
	    free(fvec);
	    exit(1);
	  }
	  d[i]=r[i][i];
	}
      }
    }
    for (i=0;i<n;i++) {
      for (sum=0.0,j=0;j<n;j++) sum += qt[i][j]*fvec[j];
      p[i]=-sum;
    }
    for (i=n-1;i>=0;i--) {
      for (sum=0.0,j=0;j<=i;j++) sum -= r[j][i]*p[j];
      g[i]=sum;
    }
    for (i=0;i<n;i++) {
      xold[i]=x[i];
      fvcold[i]=fvec[i];
    }
    fold=f;
    rsolv(r,n,d,p);
    lnsrch(n,xold,fold,g,p,x,f,stpmax,check,fvec);
    test=0.0;
    for (i=0;i<n;i++)
      if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
    if (test < TOLF) {
      check=false;
      mem->destroy_2d_double_array(qt);
      mem->destroy_2d_double_array(r);
      free(c);
      free(d);
      free(fvcold);
      free(g);
      free(p);
      free(s);
      free(t);
      free(w);
      free(xold);
      free(fvec);
      return;  
    }
    if (check) {
      if (restrt) {
	mem->destroy_2d_double_array(qt);
	mem->destroy_2d_double_array(r);
	free(c);
	free(d);
	free(fvcold);
	free(g);
	free(p);
	free(s);
	free(t);
	free(w);
	free(xold);
	free(fvec);
	return;  
      } else {
	test=0.0;
	den=MAX(f,0.5*n);
	for (i=0;i<n;i++) {
	  temp=fabs(g[i])*MAX(fabs(x[i]),1.0)/den;
	  if (temp > test) test=temp;
	}
	if (test < TOLMIN) {
	  mem->destroy_2d_double_array(qt);
	  mem->destroy_2d_double_array(r);
	  free(c);
	  free(d);
	  free(fvcold);
	  free(g);
	  free(p);
	  free(s);
	  free(t);
	  free(w);
	  free(xold);
	  free(fvec);
	  return;  
	} else restrt=true;
      }
    } else {
      restrt=false;
      test=0.0;
      for (i=0;i<n;i++) {
	temp=(fabs(x[i]-xold[i]))/MAX(fabs(x[i]),1.0);
	if (temp > test) test=temp;
      }
      if (test < TOLX)  {
	mem->destroy_2d_double_array(qt);
	mem->destroy_2d_double_array(r);
	free(c);
	free(d);
	free(fvcold);
	free(g);
	free(p);
	free(s);
	free(t);
	free(w);
	free(xold);
	free(fvec);
	return;  
      }
    }
  }
  cout<<"MAXITS exceeded in broydn"<<endl;
  mem->destroy_2d_double_array(qt);
  mem->destroy_2d_double_array(r);
  free(c);
  free(d);
  free(fvcold);
  free(g);
  free(p);
  free(s);
  free(t);
  free(w);
  free(xold);
  free(fvec);
  exit(1);
}
