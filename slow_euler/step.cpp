#include "slow_euler.hpp"
#include "slow_euler_routines.h"

void Slow_Euler::take_steps(int n){
  for(int i=0; i<n; i++){
    step();
    step_number++;
  }
}

void Slow_Euler:: go(){
  this->take_steps(inputs->steps);
}

void Slow_Euler::step(){

  FILE *f;
  double residual;
  int j, k, idx;
  int debug;

  memset(rhs, 0, dim->pts*4*sizeof(double));

  for(j=dim->nghost; j<dim->jmax+dim->nghost; j++){
  for(k=dim->nghost; k<dim->kmax+dim->nghost; k++){

    idx = j*dim->jstride + k*dim->kstride;

    slow_timestep(q[idx], grid->xy[idx], grid->xy[idx+dim->jstride], grid->xy[idx+dim->kstride], 
		  inputs->cfl, &dt[idx]);

  }
  }

  this->boundary_conditions();

  this->flux();

  residual = this->rhs_times_dt();

  if((step_number+1)%inputs->resid == 0){

    if(residual != residual){
      throw 234;
    }
    printf("%d : %12.8e\n", step_number+1, residual);
    f = fopen("residual.dat", "a");
    fprintf(f,"%d %E\n", step_number+1, residual);
    fclose(f);
  }

  if(inputs->ilhs > 0){
    this->dadi();
  }

  this->update_q();

  // if((step_number+1)%inputs->resid == 0){

  //   int idx = 190*dim->jstride + 3*dim->kstride;
  //   printf("q %e %e %e %e\n", q[idx][0], q[idx][1], q[idx][2], q[idx][3]);

  // }
}
