#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "headers.h"

void init_grid(double *x1D, double *y1D, int jtot, int ktot, int jle, 
	       double ds, double xsf, double xint){

  int j, k;
  //double xint = 1.008930411365;
  double xint2 = xint*xint;
  double xint3 = xint2*xint;
  double xint4 = xint2*xint2;
  double th = 0.12;

  double (*x)[jtot] = (double (*)[jtot])x1D;
  double (*y)[jtot] = (double (*)[jtot])y1D;

  k = 0;
  for(j=0;j<jtot;j++){

    x[k][j] = 0.5-0.5*cos(M_PI*(jle-j)/(jle));
  }
  for(j=0;j<jle;j++){
    y[k][j] = -5.0*th*(0.2969 * sqrt(xint*x[k][j]) -
		       0.1260 * xint  * x[k][j]                 -
		       0.3516 * xint2 * x[k][j]*x[k][j]         +
		       0.2843 * xint3 * x[k][j]*x[k][j]*x[k][j] -
		       0.1015 * xint4 * x[k][j]*x[k][j]*x[k][j]*x[k][j]);
  }
  for(j=jle;j<jtot;j++){
    y[k][j] = 5.0*th*(0.2969 * sqrt(xint*x[k][j]) -
		      0.1260 * xint  * x[k][j]                 -
		      0.3516 * xint2 * x[k][j]*x[k][j]         +
		      0.2843 * xint3 * x[k][j]*x[k][j]*x[k][j] -
		      0.1015 * xint4 * x[k][j]*x[k][j]*x[k][j]*x[k][j]);
  }

  //
  // Round off the trailing edge
  //
  y[0][0]      = 0.0;
  y[0][jtot-1] = 0.0;

  y[0][1]      = 0.5*y[0][1]     +0.25*(y[0][2]     +y[0][0]     );
  y[0][jtot-2] = 0.5*y[0][jtot-2]+0.25*(y[0][jtot-3]+y[0][jtot-1]);

  // // debug
  // for(j=0;j<jtot;j++){
  //   printf("%8.6lf  %8.6lf\n", x[0][j], y[0][j]);
  // }

  //
  // Go along wake direction
  //
  x[1][0] = x[0][0] + ds;
  for(k=2;k<ktot;k++){
    x[k][0] = x[k-1][0]+( x[k-1][0] - x[k-2][0] )*xsf;
  }

  double xmax = x[ktot-1][0];
  
  for(k=1;k<ktot;k++){
    y[k][0] = 0.0;
    x[k][jtot-1] = x[k][0];
    y[k][jtot-1] = y[k][0];
  }
  //
  // Fill the far-field
  //
  for(j=0;j<jtot;j++){
    x[ktot-1][j] = -xmax*cos(M_PI*(jle-j)/jle);
    y[ktot-1][j] = -xmax*sin(M_PI*(jle-j)/jle);
  }
  double dx, dy, X, Y, dr, R;
  //
  // Linear interpolation everywhere else
  //
  for(k=1;k<ktot-1;k++){
    for(j=1;j<jtot-1;j++){
      X = x[ktot-1][j]-x[0][j];
      Y = y[ktot-1][j]-y[0][j];
      R = sqrt(Y*Y+X*X);
      dr = (x[k][0]-x[k-1][0])*R/(x[ktot-1][0]-x[0][0]);
	
      dy = dr*Y/R;
      dx = dr*X/R;
      x[k][j] = x[k-1][j] + dx;
      y[k][j] = y[k-1][j] + dy;
    }
  }

  // for(j=0;j<jtot;j++){
  //   for(k=0;k<ktot;k++){
  //     printf("%8.6lf  %8.6lf\n", x[k][j], y[k][j]);
  //   }
  // }

}

void residual(double *x1D, double *y1D, double *res1D, int jtot, int ktot, int poisson){
  double (*x)[jtot]      = (double (*)[jtot])x1D;
  double (*y)[jtot]      = (double (*)[jtot])y1D;
  double (*res)[jtot][2] = (double (*)[jtot][2])res1D;

  double *psi1D = (double *)malloc(jtot*ktot*sizeof(double));
  double *phi1D = (double *)malloc(jtot*ktot*sizeof(double));

  double (*psi)[jtot] = (double (*)[jtot])psi1D;
  double (*phi)[jtot] = (double (*)[jtot])phi1D;

  double A1, A2, A3,x_eta,y_eta,x_xi,y_xi;
  double x_xi_xi, y_xi_xi, x_eta_eta, y_eta_eta;
  double tmpx, tmpy;
  int j, k;

  if(poisson){
    psi_phi(psi1D, phi1D, x1D, y1D, jtot, ktot);
  } else {
    for(k=0;k<ktot;k++){
      for(j=0;j<jtot;j++){
	phi[k][j] = 0.0;    // Homogenious Laplace Case
	psi[k][j] = 0.0;    //
      }
    }
  }

  for(k=1;k<ktot-1;k++){
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

	res[k][j][0] = ( A1*(x_xi_xi) -
			 A2*(x[k+1][j+1]-x[k+1][j-1]-x[k-1][j+1]+x[k-1][j-1])/2.0 +
			 A3*(x_eta_eta) + A1*phi[k][j]*x_xi + A3*psi[k][j]*x_eta );

	res[k][j][1] = ( A1*(y_xi_xi) -
			 A2*(y[k+1][j+1]-y[k+1][j-1]-y[k-1][j+1]+y[k-1][j-1])/2.0 +
			 A3*(y_eta_eta) + A1*phi[k][j]*y_xi + A3*psi[k][j]*y_eta );

    }
  }
}

