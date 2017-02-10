#include "adjoint.hpp"

void Adjoint::take_steps(int n){
  for(int i=0; i<n; i++){
    step();
    step_number++;
  }
}

void Adjoint:: go(){
  this->take_steps(euler->inputs->steps);
}

void Adjoint::step(){

  // FILE *f;
  // double residual;

  // memset(rhs, 0, dim->pts*4*sizeof(double));

  // //this->timestep();
  // timestep(q, grid->Sj, grid->Sk, grid->V, dim, inputs->cfl, dt);

  // this->boundary_conditions();

  // this->flux();

  // residual = this->rhs_times_dt();

  // if((step_number+1)%inputs->resid == 0){

  //   if(residual != residual){
  //     throw 234;
  //   }
  //   printf("%d : %12.8e\n", step_number+1, residual);
  //   f = fopen("residual.dat", "a");
  //   fprintf(f,"%d %E\n", step_number+1, residual);
  //   fclose(f);
  // }

  // if(inputs->ilhs > 0){
  //   this->dadi();
  // }

  // this->update_q();

  // // if((step_number+1)%inputs->resid == 0){

  // //   int idx = 190*dim->jstride + 3*dim->kstride;
  // //   printf("q %e %e %e %e\n", q[idx][0], q[idx][1], q[idx][2], q[idx][3]);

  // // }
}
