#include "adadj.hpp"
#include "adj_routines.h"

void ADadj::go(){

  int i, j, k, idx;

  // use the Q from the Euler class
  memcpy(q, euler->q, dim->pts*4*sizeof(double));

  // zero everything else
  memset(  qb, 0, 4*dim->pts*sizeof(double));
  memset(   f, 0, 4*dim->pts*sizeof(double));
  memset(  fb, 0, 4*dim->pts*sizeof(double));
  memset( rhs, 0, 4*dim->pts*sizeof(double));
  memset(rhsb, 0, 4*dim->pts*sizeof(double));
  memset( xyb, 0, 2*dim->pts*sizeof(double));
  memset(  dt, 0,   dim->pts*sizeof(double));
  memset( dtb, 0,   dim->pts*sizeof(double));

  double dummy, Jb = 1.0;

  int start, end, pj;
  start = std::max(wall->js, dim->nghost);
  end   = std::min(wall->je, dim->jmax + dim->nghost - 1);
  k = dim->nghost;
  for(j=start; j<end; j++){
    pj  = j-start;
    idx = j*dim->jstride + k*dim->kstride;

    pressure_cost_b( q[idx], qb[idx], dim, &dummy, &Jb, p_des[pj] );
    
  }

  // for(i=0; i<1; i++){
    
  //   this->step();

  // }

}

void ADadj::step(){

  double residual;

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

}
