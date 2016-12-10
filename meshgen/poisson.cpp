#include "meshgen.hpp"
#include "slor.h"
#include "ss.h"

void MeshGen::poisson(int n){

  double A1, A2, A3,x_eta,y_eta,x_xi,y_xi;
  double x_xi_xi, y_xi_xi, x_eta_eta, y_eta_eta;
  double tmpx, tmpy, res;
  int j, k;
  double w = 1.5;
  
  memset(rhs, 0, dim->pts*2*sizeof(double));
  memset(P  , 0, dim->pts*1*sizeof(double));
  memset(Q  , 0, dim->pts*1*sizeof(double));

  for(int i=0; i<n; i++){

    // middlecoff_PQ(P, Q, x, y, dim->jtot, dim->ktot);
    ss_PQ(P,Q,x,y,dim->jtot,dim->ktot,ds1,ds2);

    if((i+1)%200 == 0){
      res = residual(x,y,(double*)rhs,P,Q,dim->jtot,dim->ktot);
      printf("L2 Residual %4d : %e\n", i+1, res);
    }
    
    slor(w, x, y, (double*)rhs, P, Q, a, b, c, d, dim->jtot, dim->ktot);

  }

}
