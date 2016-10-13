#include "euler.hpp"


double Euler::rhs_times_dt(){

  int idx, j, k;

  double dt_over_volume;

  double residual = 0.0;

  for(k=dim->nghost; k<dim->kmax+dim->nghost; k++){
  for(j=dim->nghost; j<dim->jmax+dim->nghost; j++){

    idx = j*dim->jstride + k*dim->kstride;
    
    dt_over_volume = dt[idx]/grid->V[idx];

    // // if(rhs[idx][0] != rhs[idx][0]){
    // if(j==2 && k == 3){
    //   //printf("we have a NaN at (%d %d)\n", j, k);
    //   printf("%d %d: %e %e %e\n", j, k, rhs[idx][0], dt[idx], grid->V[idx]);
    // }

    rhs[idx][0] = rhs[idx][0]*dt_over_volume;
    rhs[idx][1] = rhs[idx][1]*dt_over_volume;
    rhs[idx][2] = rhs[idx][2]*dt_over_volume;
    rhs[idx][3] = rhs[idx][3]*dt_over_volume;

    residual += rhs[idx][0]*rhs[idx][0];

  }
  }

  return sqrt(residual / (dim->kmax * dim->jmax) );

}

void Euler::update_q(){

  int idx, j, k;
  int errc = 0;
  for(k=dim->nghost; k<dim->kmax+dim->nghost; k++){
  for(j=dim->nghost; j<dim->jmax+dim->nghost; j++){

    idx = j*dim->jstride + k*dim->kstride;

    q[idx][0] += rhs[idx][0];
    q[idx][1] += rhs[idx][1];
    q[idx][2] += rhs[idx][2];
    q[idx][3] += rhs[idx][3];

    // if(errc < 5 && q[idx][0] != q[idx][0]){
    //   printf("found a nan at %d %d\n", j, k);
    //   errc++;
    // }

  }
  }


}
