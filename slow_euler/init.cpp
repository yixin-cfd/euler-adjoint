#include "slow_euler.hpp"
#include "structures.h"
#include <math.h>

void Slow_Euler::init(){

  // finish setting inputs

  inputs->rho_inf = 1.0;
  inputs->u_inf   = inputs->M_inf * cos( inputs->aoa * M_PI/180.0 );
  inputs->v_inf   = inputs->M_inf * sin( inputs->aoa * M_PI/180.0 );
  inputs->p_inf   = 1.0/GAMMA;
  
  double vel_sq = inputs->M_inf * inputs->M_inf;

  inputs->e_inf   = inputs->p_inf / (GAMMA - 1.0) + 0.5 * inputs->rho_inf * vel_sq;


  printf("initializing solution\n");
  for(int i=0; i<dim->pts; i++){

    q[i][0] = inputs->rho_inf;
    q[i][1] = inputs->rho_inf * inputs->u_inf;
    q[i][2] = inputs->rho_inf * inputs->v_inf;
    q[i][3] = inputs->e_inf;

  }

}
