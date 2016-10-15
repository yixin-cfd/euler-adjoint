#include <stdio.h>
#include <math.h>

void adi(double *x1D, double *y1D, double *rhs1D, int jtot, int ktot){

  int j, k;
  double A1, A2, A3,x_eta,y_eta,x_xi,y_xi;

  double (*x)[jtot]      = (double (*)[jtot])x1D;
  double (*y)[jtot]      = (double (*)[jtot])y1D;
  double (*rhs)[jtot][2] = (double (*)[jtot][2])rhs1D;
  double a[ktot][jtot];

  double dt = 1.0; // pseudo time step

  for(j=0;j<jtot;j++){
    for(k=0;k<ktot;k++){

      // x_xi  = (x[k][j+1]-x[k][j-1])/2.0;
      // y_xi  = (y[k][j+1]-y[k][j-1])/2.0;
      x_eta = (x[k+1][j]-x[k-1][j])/2.0;
      y_eta = (y[k+1][j]-y[k-1][j])/2.0;
	
      A1 = x_eta*x_eta + y_eta*y_eta;
      // A2 = x_xi*x_eta  + y_xi*y_eta;
      // A3 = x_xi*x_xi   + y_xi*y_xi;

      a[k][j] = 1.0 + dt*A1;
    }
  }

}
