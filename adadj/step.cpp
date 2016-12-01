#include "adadj.hpp"
#include "adj_routines.h"
#include "euler_routines.h"


void ADadj::go(int nsteps){

  FILE *file;
  int i, j, k, idx;

  // use the Q from the Euler class
  memcpy(q, euler->q, dim->pts*4*sizeof(double));

  // zero everything else
  memset(rhsb, 0, 4*dim->pts*sizeof(double));
  memset( qb2, 0, 4*dim->pts*sizeof(double));
  memset(   f, 0, 4*dim->pts*sizeof(double));
  memset(  fb, 0, 4*dim->pts*sizeof(double));
  memset( rhs, 0, 4*dim->pts*sizeof(double));
  memset(  dt, 0,   dim->pts*sizeof(double));

  //
  // Derivative of cost function J wrt J is 1!
  //
  double dummy, residual, Jb = 1.0;

  //
  // Now we find derivative of cost function J wrt q
  //
  int start, end, pj;
  start = std::max(wall->js, dim->nghost);
  end   = std::min(wall->je, dim->jmax + dim->nghost - 1);
  k = dim->nghost;
  for(j=start; j<=end; j++){
    pj  = j-start;
    idx = j*dim->jstride + k*dim->kstride;

    pressure_cost_b( q[idx], qb2[idx], dim, &dummy, &Jb, p_des[pj] );
  }

  //
  // We will need the forward timestep
  //
  for(j=0; j<dim->jtot; j++){
  for(k=0; k<dim->ktot; k++){
    idx = j*dim->jstride + k*dim->kstride;

    ad_timestep(q[idx], grid->xy[idx], grid->xy[idx+dim->jstride], grid->xy[idx+dim->kstride], 
		euler->inputs->cfl, &dt[idx]);

  }
  }

  printf("nsteps = %d\n", nsteps);
  //
  // Main iteration Loop
  //
  file = fopen("res_adj.dat", "w");
  fclose(file);
    
  for(i=0; i<nsteps; i++){
    residual = this->step();
  
    step_number++;

    if(step_number % euler->inputs->resid == 0){
      residual = sqrt(residual/(dim->jtot*dim->ktot));
      printf("%4d : %15.6e\n", step_number, residual);
      file = fopen("res_adj.dat", "a");
      fprintf(file,"%d %E\n", step_number, residual);
      fclose(file);
    }

  }


}

double ADadj::step(){
  
  double residual=0.0;
  
  int j, k, idx;
  int jstride = dim->jstride;
  int kstride = dim->kstride;
  
  memset(  qb, 0, 4*dim->pts*sizeof(double));
  memset( dtb, 0,   dim->pts*sizeof(double));
    
  this->flux(false);
  this->boundary_conditions(false);

  for(j=0; j<dim->jtot-1; j++){
  for(k=0; k<dim->ktot-1; k++){

    idx = j*dim->jstride + k*dim->kstride;
    
    ad_timestep_b(q[idx], qb[idx], grid->xy[idx], grid->xy[idx+jstride], grid->xy[idx+kstride],
    		  euler->inputs->cfl, &dt[idx], &dtb[idx]);

  }
  }

  // this->boundary_conditions(false);

  // j = 43; k = 66;
  // idx = j*dim->jstride + k*dim->kstride;
  // printf("__ %15.8e %15.8e %15.8e %15.8e\n",
  // 	 rhsb[idx][0],rhsb[idx][1],rhsb[idx][2],rhsb[idx][3]);
  // printf("__ %15.8e %15.8e %15.8e %15.8e\n",
  // 	 qb2[idx][0],qb2[idx][1],qb2[idx][2],qb2[idx][3]);

  // for(j=dim->nghost; j<dim->jmax+dim->nghost; j++){
  // for(k=dim->nghost; k<dim->kmax+dim->nghost; k++){
  for(j=0; j<dim->jtot; j++){
  for(k=0; k<dim->ktot; k++){

    idx = j*dim->jstride + k*dim->kstride;

    // add contribution from cost function
    qb[idx][0] = -qb2[idx][0] + qb[idx][0];
    qb[idx][1] = -qb2[idx][1] + qb[idx][1];
    qb[idx][2] = -qb2[idx][2] + qb[idx][2];
    qb[idx][3] = -qb2[idx][3] + qb[idx][3];

    // update the residual
    rhsb[idx][0] = rhsb[idx][0] + qb[idx][0]*dt[idx];
    rhsb[idx][1] = rhsb[idx][1] + qb[idx][1]*dt[idx];
    rhsb[idx][2] = rhsb[idx][2] + qb[idx][2]*dt[idx];
    rhsb[idx][3] = rhsb[idx][3] + qb[idx][3]*dt[idx];

    residual += qb[idx][0]*qb[idx][0];
    residual += qb[idx][1]*qb[idx][1];
    residual += qb[idx][2]*qb[idx][2];
    residual += qb[idx][3]*qb[idx][3];

    if(residual != residual){
      printf("nan at %d %d\n", j, k);
      throw 2341;
    }

  }
  }

  return residual;
  
}

void ADadj::check(){

  int j, k, idx;
  int jstride = dim->jstride;
  int kstride = dim->kstride;
  
  memset( xyb, 0, 2*dim->pts*sizeof(double));

  this->boundary_conditions(true);

  this->flux(true);

  // for(j=dim->nghost; j<dim->jmax+dim->nghost; j++){
  // for(k=dim->nghost; k<dim->kmax+dim->nghost; k++){

  //   idx = j*dim->jstride + k*dim->kstride;

  //   ad_timestep_bx(q[idx], qb[idx], grid->xy[idx], xyb[idx],
  //   		  grid->xy[idx+jstride], xyb[idx+jstride],
  //   		  grid->xy[idx+kstride], xyb[idx+kstride],
  //   		  euler->inputs->cfl, &dt[idx], &dtb[idx]);
  // }
  // }
  
}
