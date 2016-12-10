#include "meshgen.hpp"

#define XSF 1.12

void MeshGen::init(double dist){

  int jstride = dim->jstride;
  int kstride = dim->kstride;
  int jtot = dim->jtot;
  int ktot = dim->ktot;

  double ds;
  if(dist > 0.0 && dist < 0.1){
    ds = dist;
  } else {
    ds = 0.5*(x[0] - x[2*jstride]);
  }

  int idx, j, k, jle;

  jle = (dim->jtot - 1)/2;

  for(j=0; j<dim->pts; j++){
    x[j] = x[j] - 0.5;
  }
  
  //
  // Go along wake direction
  //
  j = 0;
  x[kstride] = x[0] + ds;
  for(k=2;k<ktot;k++){
    idx = j*jstride + k*kstride;
    x[idx] = x[idx - kstride]+( x[idx-kstride] - x[idx - 2*kstride] )*XSF;
  }

  double xmax = x[(ktot-1)*kstride];

  // printf("warning: not using dist. How close? %f vs %f\n", dist, xmax);

  j = jtot-1;
  for(k=1;k<ktot;k++){
    idx = j*jstride + k*kstride;
    y[k*kstride] = 0.0;
    x[idx] = x[k*kstride];
    y[idx] = y[k*kstride];
  }

  //
  // Fill the far-field
  //
  k = ktot-1;
  for(j=0;j<jtot;j++){
    idx = j*jstride + k*kstride;
    x[idx] = -xmax*cos(M_PI*(jle-j)/jle);
    y[idx] = -xmax*sin(M_PI*(jle-j)/jle);
  }

  double dx, dy, X, Y, dr, R;
  //
  // Linear interpolation everywhere else
  //
  for(j=1;j<jtot-1;j++){

    X = x[j*jstride + (ktot-1)*kstride]-x[j*jstride];
    Y = y[j*jstride + (ktot-1)*kstride]-y[j*jstride];
    R = sqrt(Y*Y+X*X);
    
    for(k=1;k<ktot-1;k++){

      idx = j*jstride + k*kstride;
      dr = R*(x[k*kstride]-x[(k-1)*kstride])/(x[(ktot-1)*kstride]-x[0]);

      dy = dr*Y/R;
      dx = dr*X/R;
      x[idx] = x[idx - kstride] + dx;
      y[idx] = y[idx - kstride] + dy;
    }
  }

  for(j=0; j<dim->pts; j++){
    x[j] = x[j] + 0.5;
  }


}
