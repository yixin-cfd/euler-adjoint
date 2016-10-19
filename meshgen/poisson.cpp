#include "meshgen.hpp"
#include "slor.h"


void MeshGen::poisson(int n){

  double A1, A2, A3,x_eta,y_eta,x_xi,y_xi;
  double x_xi_xi, y_xi_xi, x_eta_eta, y_eta_eta;
  double tmpx, tmpy, res;
  int j, k;
  double w = 1.6;
  
  memset(rhs, 0, dim->pts*2*sizeof(double));
  memset(P  , 0, dim->pts*1*sizeof(double));
  memset(Q  , 0, dim->pts*1*sizeof(double));

  for(int i=0; i<n; i++){

    fill_PQ(P, Q, x, y, dim->jtot, dim->ktot);

    if(i%10 == 1){
      res = residual(x,y,(double*)rhs,P,Q,dim->jtot,dim->ktot);
      printf("residual : %e\n", res);
    }
    
    slor(w, x, y, (double*)rhs, P, Q, a, b, c, d, dim->jtot, dim->ktot);

  }

}
