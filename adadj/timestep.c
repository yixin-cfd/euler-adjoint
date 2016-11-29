#include "structures.h"
#include <math.h>
#include <stdlib.h>

void timestep(double q[4], double xy1[2], double xy2[2],  double xy3[2],
	      double cfl, double *dt){

  // xy1  is j  , k   pt
  // xy2  is j+1, k   pt
  // xy3  is j  , k+1 pt
  double Sk[2], Sj[2];
  
  Sk[0] = - xy2[1] + xy1[1];
  Sk[1] =   xy2[0] - xy1[0];
  Sj[0] =   xy3[1] - xy1[1];
  Sj[1] = - xy3[0] + xy1[0];

  //
  // Constant-CFL timestep
  //
  int idx, j, k;
  double u, v, p, c2, irho;
  double uu, vv;
  double xs2, ys2, xsc, ysc, eigmax;

  irho = 1.0 / q[0];

  u  = q[1]*irho;
  v  = q[2]*irho;
  p  = (GAMMA-1.0)*(q[3] - 0.5*q[0]*(u*u + v*v));
  c2 = GAMMA*p*irho;

  // contravariant velocities in coordinate directions
  uu = u*Sj[0] + v*Sj[1];
  vv = u*Sk[0] + v*Sk[1];

  // face sizes squared
  xs2 = (Sj[0]*Sj[0] + Sj[1]*Sj[1]);
  ys2 = (Sk[0]*Sk[0] + Sk[1]*Sk[1]);
    
  xsc = sqrt(c2*xs2);
  ysc = sqrt(c2*ys2);

  eigmax = abs(uu) + xsc + abs(vv) + ysc;
  //eigmax = fmax( abs(uu) + xsc, abs(vv) + ysc );

  // dt[0] = cfl * V[0] / eigmax;

  // actually dt over volume since later we would divide dt by volume anyway
  dt[0] = cfl / eigmax;

}
