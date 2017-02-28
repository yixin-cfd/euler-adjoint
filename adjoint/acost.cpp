#include "adjoint.hpp"

// Cost function contribution to
void Adjoint::cost(double (*rhs)[4]){

  int j, k, idx, start, end;
  double p, pd, u, v;
  double gm1 = GAMMA-1.0;

  k = dim->nghost;

  start = std::max(wall->js, dim->nghost);
  end   = std::min(wall->je, dim->jmax + dim->nghost - 1);

  for(j=start; j<=end; j++){

    idx = j*dim->jstride + k*dim->kstride;
    u   = q[idx][1]/q[idx][0];
    v   = q[idx][2]/q[idx][0];
    p   = (GAMMA - 1.0)*(q[idx][3] - 0.5*q[idx][0]*(u*u + v*v));
    pd  = this->p_des[j-start];

    rhs[idx][0] = 0.5*(p-pd)*(gm1)*(u*u + v*v);
    rhs[idx][1] =    -(p-pd)*(gm1)*(u);
    rhs[idx][2] =    -(p-pd)*(gm1)*(v);
    rhs[idx][3] =     (p-pd)*(gm1);
    
  }

}
