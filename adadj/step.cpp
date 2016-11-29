#include "adadj.hpp"
#include "adj_routines.h"
#include "euler_routines.h"


void ADadj::go(){

  int i, j, k, idx;

  // use the Q from the Euler class
  memcpy(q, euler->q, dim->pts*4*sizeof(double));

  // zero everything else
  memset( qb2, 0, 4*dim->pts*sizeof(double));
  memset(   f, 0, 4*dim->pts*sizeof(double));
  memset(  fb, 0, 4*dim->pts*sizeof(double));
  memset( rhs, 0, 4*dim->pts*sizeof(double));
  memset( xyb, 0, 2*dim->pts*sizeof(double));
  memset(  dt, 0,   dim->pts*sizeof(double));

  //
  // Derivative of cost function J wrt J is 1!
  //
  double dummy, Jb = 1.0;

  //
  // Now we find derivative of cost function J wrt q
  //
  int start, end, pj;
  start = std::max(wall->js, dim->nghost);
  end   = std::min(wall->je, dim->jmax + dim->nghost - 1);
  k = dim->nghost;
  for(j=start; j<end; j++){
    pj  = j-start;
    idx = j*dim->jstride + k*dim->kstride;

    pressure_cost_b( q[idx], qb2[idx], dim, &dummy, &Jb, p_des[pj] );
  }

  //
  // We will need the forward timestep
  //
  timestep(q, grid->Sj, grid->Sk, grid->V, dim, euler->inputs->cfl, dt);

  //
  // Main iteration Loop
  //
  for(i=0; i<1; i++){
    this->step();
  }
}

void ADadj::step(){

  double residual=0.0;

  int i, j, k;

  memset(rhsb, 0, 4*dim->pts*sizeof(double));
  memset(  qb, 0, 4*dim->pts*sizeof(double));
  memset( dtb, 0,   dim->pts*sizeof(double));

  this->boudary_conditions();

  this->flux();

  for(j=dim->nghost; j<dim->jmax+dim->nghost; j++){
  for(k=dim->nghost; k<dim->kmax+dim->nghost; k++){

    idx = j*dim->jstride + k*dim->kstride;

    // timestep(q[idx], grid->xy[idx], grid->xy[idx+jstride], grid->xy[idx+kstride],
    // 	     dim, euler->inputs->cfl, dt);

    timestep_b(q[idx], qb[idx], grid->xy[idx], xyb[idx],
	       grid->xy[idx+jstride], xyb[idx+jstride],
	       grid->xy[idx+kstride], xyb[idx+kstride],
	       dim, euler->inputs->cfl, &dt[idx], &dtb[idx]);

    // add contribution from
    qb[idx][0] = qb[idx][0] - qb2[idx][0];
    qb[idx][1] = qb[idx][1] - qb2[idx][1];
    qb[idx][2] = qb[idx][2] - qb2[idx][2];
    qb[idx][3] = qb[idx][3] - qb2[idx][3];
    
    // update the residual
    resb[idx][0] = resb[idx][0] + qb[idx][0]*dt[idx];
    resb[idx][1] = resb[idx][1] + qb[idx][1]*dt[idx];
    resb[idx][2] = resb[idx][2] + qb[idx][2]*dt[idx];
    resb[idx][3] = resb[idx][3] + qb[idx][3]*dt[idx];
    
    residual += qb[idx][0]*qb[idx][0];

  }
  }

  residual = sqrt(residual/(dim->jtot*dim->ktot));
  

  
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
