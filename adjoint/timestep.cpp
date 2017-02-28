#include "adjoint.hpp"


void Adjoint::timestep(){

  int idx, j, k;
  double u, v, p, c2, irho;
  double uu, vv;
  double xs2, ys2, xsc, ysc, eigmax;

  for(j=0; j<dim->jtot; j++){
  for(k=0; k<dim->ktot; k++){

    idx = j*dim->jstride + k*dim->kstride;

    irho = 1.0 / q[idx][0];

    u  = q[idx][1]*irho;
    v  = q[idx][2]*irho;
    p  = (GAMMA-1.0)*(q[idx][3] - 0.5*q[idx][0]*(u*u + v*v));
    c2 = GAMMA*p*irho;

    // contravariant velocities in coordinate directions
    uu = u*grid->Sj[idx][0] + v*grid->Sj[idx][1];
    vv = u*grid->Sk[idx][0] + v*grid->Sk[idx][1];

    // face sizes squared
    xs2 = (grid->Sj[idx][0]*grid->Sj[idx][0] + grid->Sj[idx][1]*grid->Sj[idx][1]);
    ys2 = (grid->Sk[idx][0]*grid->Sk[idx][0] + grid->Sk[idx][1]*grid->Sk[idx][1]);
    
    xsc = sqrt(c2*xs2);
    ysc = sqrt(c2*ys2);

    eigmax = std::abs(uu) + xsc + std::abs(vv) + ysc;
    //eigmax = fmax( abs(uu) + xsc, abs(vv) + ysc );

    // dt[0] = cfl * V[0] / eigmax;

    // actually dt over volume since later we would divide dt by volume anyway
    dt[idx] = euler->inputs->cfl / eigmax;

  }
  }


}
