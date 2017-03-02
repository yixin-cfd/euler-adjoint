#include "adjoint.hpp"

void Adjoint::take_steps(int nsteps){

  FILE *file;
  double residual, liftd;

  file = fopen("res_adj.dat", "w");
  fclose(file);

  cost(rhs0);

  timestep();

  for(int i=0; i<nsteps; i++){
    residual = step();
    step_number++;

    // if(step_number % euler->inputs->resid == 0 || step_number == nsteps){
    //   residual = sqrt(residual/(dim->jtot*dim->ktot));
    //   liftd    = 1.0;
    //   printf("%5d  %15.6e  %15.8e\n", step_number, residual, liftd);
    //   file = fopen("res_adj.dat", "a");
    //   fprintf(file, "%5d  %15.6e  %15.8e\n", step_number, residual, liftd);
    //   fclose(file);
    //   if(residual > 10) return;
    // }

    printf("Residual is %25.16e\n", residual);
  }
}

void Adjoint:: go(){
  this->take_steps(euler->inputs->steps);
}

double Adjoint::step(){

  // FILE *f;
  int j, k, idx;
  double residual = 0.0;

  memset( rhs, 0, 4*dim->pts*sizeof(double));

  this->aflux();


  // j = 5; k = 1;
  // idx = j*dim->jstride + k*dim->kstride;
  // printf("__ %20.16e %20.16e %20.16e %20.16e \n", 
  // 	 rhs[idx][0], rhs[idx][1], rhs[idx][2], rhs[idx][3]);
  // printf("__ %20.16e %20.16e %20.16e %20.16e \n", 
  // 	     psi[idx][0], psi[idx][1], psi[idx][2], psi[idx][3]);

  this->boundary_conditions();

  for(j=0; j<dim->jtot; j++){
  for(k=0; k<dim->ktot; k++){

    idx = j*dim->jstride + k*dim->kstride;

    // if(j == 5 && k == 1){
    //   printf("__ %20.16e %20.16e %20.16e %20.16e \n", 
    //   	     rhs[idx][0], rhs[idx][1], rhs[idx][2], rhs[idx][3]);
    //   // printf("__ %20.16e %20.16e %20.16e %20.16e \n", 
    //   // 	     psi[idx][0], psi[idx][1], psi[idx][2], psi[idx][3]);
    // }

    // add contribution from cost function
    rhs[idx][0] = rhs0[idx][0] + rhs[idx][0];
    rhs[idx][1] = rhs0[idx][1] + rhs[idx][1];
    rhs[idx][2] = rhs0[idx][2] + rhs[idx][2];
    rhs[idx][3] = rhs0[idx][3] + rhs[idx][3];

    psi[idx][0] = psi[idx][0] - rhs[idx][0] * dt[idx];
    psi[idx][1] = psi[idx][1] - rhs[idx][1] * dt[idx];
    psi[idx][2] = psi[idx][2] - rhs[idx][2] * dt[idx];
    psi[idx][3] = psi[idx][3] - rhs[idx][3] * dt[idx];

    residual += rhs[idx][0]*rhs[idx][0];
    residual += rhs[idx][1]*rhs[idx][1];
    residual += rhs[idx][2]*rhs[idx][2];
    residual += rhs[idx][3]*rhs[idx][3];

    if(residual != residual){
      printf("nan at %d %d\n", j, k);
      throw 2341;
    }

  }
  }

  return residual;

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
