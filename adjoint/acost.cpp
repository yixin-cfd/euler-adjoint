#include "adjoint.hpp"

// Cost function contribution to
void Adjoint::cost(double (*rhs)[4]){

  int i, j, k, idx, start, end;
  double p, pd;

  k = dim->nghost;

  start = std::max(wall->js, dim->nghost);
  end   = std::min(wall->je, dim->jmax + dim->nghost - 1);


    for(j=bc.js; j<=bc.je; j++){

      idx = j*dim->jstride + k*dim->kstride;
      p   = (GAMMA - 1.0)*(q[idx][3] - 0.5*q[idx][0]*(u*u + v*v));
      pd  = this->pdesired[j]

      rhs[idx][0] = 


    }


  double p;
  double dynp = 0.5*inputs->rho_inf*inputs->M_inf*inputs->M_inf;
  double tempb;
  double tempb0;
  // pressure at the wall is approximately the pressure at the cell right above the wall
  p = ((1.4-1.0)*(q[3]-0.5*(q[1]*q[1]+q[2]*q[2])/q[0])-inputs->p_inf)/dynp;
  pb = 2*(p-p_desired);
  tempb = 0.4*pb/dynp;
  tempb0 = -(0.5*tempb/q[0]);
  qb[3] = qb[3] + tempb;
  qb[1] = qb[1] + 2*q[1]*tempb0;
  qb[2] = qb[2] + 2*q[2]*tempb0;
  qb[0] = qb[0] - (q[1]*q[1]+q[2]*q[2])*tempb0/q[0];


}
