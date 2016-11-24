#include "structures.h"
#include "euler_routines.h"
#include <math.h>
#include <stdlib.h>

void timestep(double (*q)[4], double (*Sj)[2], double (*Sk)[2],  
	      double *V, Dim *dim, double cfl, double *dt){
  //
  // Constant-CFL timestep
  //
  int idx, j, k;
  double u, v, p, c2, irho;
  double uu, vv;
  double xs2, ys2, xsc, ysc, eigmax;

  for(k=dim->nghost; k<dim->kmax+dim->nghost; k++){
  for(j=dim->nghost; j<dim->jmax+dim->nghost; j++){

    idx = j*dim->jstride + k*dim->kstride;

    irho = 1.0 / q[idx][0];

    u  = q[idx][1]*irho;
    v  = q[idx][2]*irho;
    p  = (GAMMA-1.0)*(q[idx][3] - 0.5*q[idx][0]*(u*u + v*v));
    c2 = GAMMA*p*irho;

    // contravariant velocities in coordinate directions
    uu = u*Sj[idx][0] + v*Sj[idx][1];
    vv = u*Sk[idx][0] + v*Sk[idx][1];

    // face sizes squared
    xs2 = (Sj[idx][0]*Sj[idx][0] + Sj[idx][1]*Sj[idx][1]);
    ys2 = (Sk[idx][0]*Sk[idx][0] + Sk[idx][1]*Sk[idx][1]);
    
    xsc = sqrt(c2*xs2);
    ysc = sqrt(c2*ys2);

    //eigmax = std::abs(uu) + xsc + std::abs(vv) + ysc;
    eigmax = fmax( abs(uu) + xsc, abs(vv) + ysc );

    dt[idx] = cfl * V[idx] / eigmax;

  }
  }

}
