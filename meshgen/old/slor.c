#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "headers.h"

void psi_phi(double *psi1D, double *phi1D, double *x1D, double *y1D, int jtot, int ktot){

  double x_xi, y_xi, x_xi_xi, y_xi_xi;
  double x_eta, y_eta, x_eta_eta, y_eta_eta;

  int j, k;

  double (*x)[jtot]      = (double (*)[jtot])x1D;
  double (*y)[jtot]      = (double (*)[jtot])y1D;

  double (*phi)[jtot] = (double (*)[jtot])phi1D;
  double (*psi)[jtot] = (double (*)[jtot])psi1D;

  double jfac, kfac;

  //
  // J=0, JMAX boundary psi (xi constant)
  //
  for(k=1;k<ktot-1;k++){
    for(j=0;j<jtot;j+=jtot-1){
      x_eta     = (x[k+1][j]-x[k-1][j])/2.0;
      y_eta     = (y[k+1][j]-y[k-1][j])/2.0;
	
      x_eta_eta = x[k+1][j]-2.0*x[k][j]+x[k-1][j];
      y_eta_eta = y[k+1][j]-2.0*y[k][j]+y[k-1][j];

      if(fabs(y_eta) > fabs(x_eta))
	psi[k][j] = -y_eta_eta/y_eta;
      else
	psi[k][j] = -x_eta_eta/x_eta;

      // printf("psi %3d %3d: %16.8e\n", j, k, psi[k][j]);
    }
  }
  //
  // K=0, KMAX boundary phi (eta constant)
  //
  for(j=1;j<jtot-1;j++){
    for(k=0;k<ktot;k+=ktot-1){
      x_xi      = (x[k][j+1]-x[k][j-1])/2.0;
      y_xi      = (y[k][j+1]-y[k][j-1])/2.0;

      x_xi_xi   = x[k][j-1]-2.0*x[k][j]+x[k][j+1];
      y_xi_xi   = y[k][j-1]-2.0*y[k][j]+y[k][j+1];

      //if(y_xi*y_xi > x_xi*x_xi){
      if(fabs(y_xi) > fabs(x_xi)){
	//printf("^");
	phi[k][j] = -y_xi_xi/y_xi;
      } else {
	//printf("*");
	phi[k][j] = -x_xi_xi/x_xi;
      }

      // printf("phi %3d %3d: %16.8e\t", j, k, phi[k][j]);
      // printf("%lf %lf", x_xi, y_xi);
      // printf("\n");
      
    }
  }
  //
  // Now interpolate
  //
  for(k=1;k<ktot-1;k++){
    for(j=1;j<jtot-1;j++){

      jfac = double(j)/(jtot-1);
      kfac = double(k)/(ktot-1);
      
      psi[k][j] = psi[k][0] + jfac*(psi[k][jtot-1] - psi[k][0]);
      phi[k][j] = phi[0][j] + kfac*(phi[ktot-1][j] - phi[0][j]);
      // psi[k][j] = 0.0;
      // phi[k][j] = 0.0;
    }
  }
}

void slor(double w, double *x1D, double *y1D, int jtot, int ktot, int poisson){

  int j, k, dir;
  double A1, A2, A3,x_eta,y_eta,x_xi,y_xi;
  double y_eta_eta, x_eta_eta, y_xi_xi, x_xi_xi;

  double (*x)[jtot]      = (double (*)[jtot])x1D;
  double (*y)[jtot]      = (double (*)[jtot])y1D;

  double *a = (double *)malloc(jtot*sizeof(double));
  double *b = (double *)malloc(jtot*sizeof(double));
  double *c = (double *)malloc(jtot*sizeof(double));
  double *d = (double *)malloc(jtot*sizeof(double));

  double *psi1D = (double *)malloc(jtot*ktot*sizeof(double));
  double *phi1D = (double *)malloc(jtot*ktot*sizeof(double));

  double (*psi)[jtot] = (double (*)[jtot])psi1D;
  double (*phi)[jtot] = (double (*)[jtot])phi1D;
  double res = 0.0;

  double J;

  if(poisson){
    psi_phi(psi1D, phi1D, x1D, y1D, jtot, ktot);
  } else {
    for(k=0;k<ktot;k++){
      for(j=0;j<jtot;j++){  //
	phi[k][j] = 0.0;    // Homogenious Laplace Case
	psi[k][j] = 0.0;    //
      }
    }
  }
  
  for(k=1;k<ktot-1;k++){
    // loop through x and y directions
    for(dir=0; dir<2; dir++){
      //
      // Boundary Conditions
      //
      a[0] = 0;
      b[0] = 1;
      c[0] = 0;
      d[0] = 0;

      for(j=1;j<jtot-1;j++){

	x_xi      = (x[k][j+1]-x[k][j-1])/2.0;
	y_xi      = (y[k][j+1]-y[k][j-1])/2.0;

	x_eta     = (x[k+1][j]-x[k-1][j])/2.0;
	y_eta     = (y[k+1][j]-y[k-1][j])/2.0;
	
	x_xi_xi   = x[k][j-1]-2.0*x[k][j]+x[k][j+1];
	y_xi_xi   = y[k][j-1]-2.0*y[k][j]+y[k][j+1];

	x_eta_eta = x[k+1][j]-2.0*x[k][j]+x[k-1][j];
	y_eta_eta = y[k+1][j]-2.0*y[k][j]+y[k-1][j];
	
	A1 = x_eta*x_eta + y_eta*y_eta;
	A2 = x_xi*x_eta  + y_xi*y_eta;
	A3 = x_xi*x_xi   + y_xi*y_xi;

	a[j] = A1;
	b[j] = -2.0*A1-2.0*A3;
	c[j] = A1;

	if(dir == 0){

	  d[j] = -w*( A1*(x_xi_xi) -
		      A2*0.5*(x[k+1][j+1]-x[k+1][j-1]-x[k-1][j+1]+x[k-1][j-1]) +
		      A3*(x_eta_eta) + A1*phi[k][j]*x_xi + A3*psi[k][j]*x_eta );

	} else {
	  d[j] = -w*( A1*(y_xi_xi) -
		      A2*0.5*(y[k+1][j+1]-y[k+1][j-1]-y[k-1][j+1]+y[k-1][j-1]) +
		      A3*(y_eta_eta) + A1*phi[k][j]*y_xi + A3*psi[k][j]*y_eta );

	}

	//res += d[j]*d[j];
      }

      a[jtot-1] = 0;
      b[jtot-1] = 1;
      c[jtot-1] = 0;
      d[jtot-1] = 0;

      tri(a,b,c,d,jtot);
      
      if(dir == 0){ // x-dir
	for(j=1;j<jtot-1;j++){
	  x[k][j] += d[j];
	}
      } else {
	for(j=1;j<jtot-1;j++){
	  y[k][j] += d[j];
	}
      }
    }
  }

  //printf("Norm Residual this round: %16.8e\n", sqrt(res/(jtot*ktot)));

  free(a);
  free(b);
  free(c);
  free(d);
  free(phi1D);
  free(psi1D);
}
