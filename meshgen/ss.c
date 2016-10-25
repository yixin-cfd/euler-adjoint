#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define AAA 0.45
#define BBB 0.45
#define CCC 0.45
#define DDD 0.45

#define OMEGA 0.1

void ss_PQ(double *P1D, double *Q1D, double *x1D, double *y1D, 
	   int jtot, int ktot){

  double tmp, alpha, beta, gamma, iJ,x_xi,y_xi;
  double x_xi_xi, y_xi_xi, x_eta_eta, y_eta_eta;
  double x_eta_xi, y_eta_xi;
  double sn, R1, R2, R3, R4, newP, newQ;

  int j, k;

  double (*x)[jtot]      = (double (*)[jtot])x1D;
  double (*y)[jtot]      = (double (*)[jtot])y1D;

  double (*P)[jtot] = (double (*)[jtot])P1D;
  double (*Q)[jtot] = (double (*)[jtot])Q1D;

  double *p, *q, *r, *s;
  double *x_eta, *y_eta;

  p     = (double*)malloc(jtot*sizeof(double));
  q     = (double*)malloc(jtot*sizeof(double));
  r     = (double*)malloc(jtot*sizeof(double));
  s     = (double*)malloc(jtot*sizeof(double));
  x_eta = (double*)malloc(jtot*sizeof(double));
  y_eta = (double*)malloc(jtot*sizeof(double));
  
  // boundaries
  j = 0;
  x_eta[j] = x[1][j]-x[0][j];
  y_eta[j] = y[1][j]-y[0][j];
  j = jtot-1;
  x_eta[j] = x[1][j]-x[0][j];
  y_eta[j] = y[1][j]-y[0][j];

  // fill the rest of the eta derivatives
  k = 0;
  for(j=1;j<jtot-1;j++){
    sn       = x[1][0]-x[0][0]; // because we forced this point along
				// the wake
    x_xi     = (x[k][j+1]-x[k][j-1])/2.0;
    y_xi     = (y[k][j+1]-y[k][j-1])/2.0;
    gamma    = x_xi*x_xi   + y_xi*y_xi;
    tmp      = 1.0/sqrt(gamma);
    x_eta[j] = sn * (-y_xi) * tmp;
    y_eta[j] = sn * ( x_xi) * tmp;
  }
  
  k = 0;
  for(j=1;j<jtot-1;j++){

    x_xi      = 0.5*(x[k][j+1]-x[k][j-1]);
    y_xi      = 0.5*(y[k][j+1]-y[k][j-1]);

    alpha = x_eta[j]*x_eta[j] + y_eta[j]*y_eta[j];
    beta  = x_xi*x_eta[j]  + y_xi*y_eta[j];
    gamma = x_xi*x_xi   + y_xi*y_xi;
    iJ    =   1.0/(x_xi *y_eta[j] - y_xi*x_eta[j] );

    x_xi_xi   = x[k][j-1]-2.0*x[k][j]+x[k][j+1];
    y_xi_xi   = y[k][j-1]-2.0*y[k][j]+y[k][j+1];
    
    x_eta_eta = .5*(-7.0*x[0][j] + 8.0*x[1][j] - x[2][j]) - 3.0*x_eta[j];
    y_eta_eta = .5*(-7.0*y[0][j] + 8.0*y[1][j] - y[2][j]) - 3.0*y_eta[j];
    x_eta_xi  = .5*(x_eta[j+1]-x_eta[j-1]);
    y_eta_xi  = .5*(y_eta[j+1]-y_eta[j-1]);

    R1 = -(alpha*x_xi_xi -2.0*beta*x_eta_xi + gamma*x_eta_eta)*iJ*iJ;
    R2 = -(alpha*y_xi_xi -2.0*beta*y_eta_xi + gamma*y_eta_eta)*iJ*iJ;

    p[j] = ( y_eta[j]*R1 - x_eta[j]*R2)*iJ;
    q[j] = (-y_xi*R1  + x_xi*R2)*iJ;
  }

  // done with eta-min boundary (k=0)

  k = ktot-1;
  // j-boundaries
  j = 0;
  x_eta[j] = x[k][j]-x[k-1][j];
  y_eta[j] = y[k][j]-y[k-1][j];
  j = jtot-1;
  x_eta[j] = x[k][j]-x[k-1][j];
  y_eta[j] = y[k][j]-y[k-1][j];

  // fill the rest of the eta derivatives
  for(j=1;j<jtot-1;j++){
    sn       = x[k][0]-x[k-1][0]; // because we forced this point
				  // along the wake
    x_xi     = (x[k][j+1]-x[k][j-1])/2.0;
    y_xi     = (y[k][j+1]-y[k][j-1])/2.0;
    gamma    = x_xi*x_xi   + y_xi*y_xi;
    tmp      = 1.0/sqrt(gamma);
    x_eta[j] = sn * (-y_xi) * tmp;
    y_eta[j] = sn * ( x_xi) * tmp;
  }

  for(j=1;j<jtot-1;j++){

    x_xi      = 0.5*(x[k][j+1]-x[k][j-1]);
    y_xi      = 0.5*(y[k][j+1]-y[k][j-1]);

    alpha = x_eta[j]*x_eta[j] + y_eta[j]*y_eta[j];
    beta  = x_xi*x_eta[j]  + y_xi*y_eta[j];
    gamma = x_xi*x_xi   + y_xi*y_xi;
    iJ    =   1.0/(x_xi *y_eta[j] - y_xi*x_eta[j] );
    
    x_xi_xi   = x[k][j-1]-2.0*x[k][j]+x[k][j+1];
    y_xi_xi   = y[k][j-1]-2.0*y[k][j]+y[k][j+1];
    
    x_eta_eta = .5*(-7.0*x[k][j] + 
		    8.0*x[k-1][j] - x[k-2][j]) + 3.0*x_eta[j];
    y_eta_eta = .5*(-7.0*y[k][j] + 
		    8.0*y[k-1][j] - y[k-2][j]) + 3.0*y_eta[j];

    x_eta_xi = 0.5*(x_eta[j+1]-x_eta[j-1]);
    y_eta_xi = 0.5*(y_eta[j+1]-y_eta[j-1]);

    R3 = -(alpha*x_xi_xi -2.0*beta*x_eta_xi + gamma*x_eta_eta)*iJ*iJ;
    R4 = -(alpha*y_xi_xi -2.0*beta*y_eta_xi + gamma*y_eta_eta)*iJ*iJ;

    r[j] = ( y_eta[j]*R3 - x_eta[j]*R4)*iJ;
    s[j] = (-y_xi*R3  + x_xi*R4)*iJ;
  }

  // ok now we have p,q,r,s for all j
  for(k=1; k<ktot-1; k++){
  for(j=1; j<jtot-1; j++){

    newP = p[j]*exp(-AAA*k) + r[j]*exp(-CCC*(ktot-1-k));
    newQ = q[j]*exp(-BBB*k) + s[j]*exp(-DDD*(ktot-1-k));

    P[k][j] = P[k][j]*(1-OMEGA) + OMEGA*newP;
    Q[k][j] = Q[k][j]*(1-OMEGA) + OMEGA*newQ;

  }
  }


  free(p);
  free(q);

}
