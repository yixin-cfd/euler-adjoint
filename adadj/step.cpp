#include "adadj.hpp"
#include "adj_routines.h"
#include "euler_routines.h"


void ADadj::go(int nsteps){

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
  for(j=0; j<dim->jtot; j++){
  for(k=0; k<dim->ktot; k++){
    idx = j*dim->jstride + k*dim->kstride;

    ad_timestep(q[idx], grid->xy[idx], grid->xy[idx+dim->jstride], grid->xy[idx+dim->kstride], 
		euler->inputs->cfl, &dt[idx]);

    if(j==8 && k ==1){
      printf("timestep is: %f\n", dt[idx]*grid->V[idx]);
    }
  }
  }

  //
  // Main iteration Loop
  //
  for(i=0; i<nsteps; i++){
    this->step();
  }
}

void ADadj::step(){
  
  double residual=0.0;
  
  int j, k, idx;
  int jstride = dim->jstride;
  int kstride = dim->kstride;
  
  memset(  qb, 0, 4*dim->pts*sizeof(double));
  memset( dtb, 0,   dim->pts*sizeof(double));
    
  this->flux(false);
  this->boundary_conditions(false);

  for(j=0; j<dim->jtot; j++){
  for(k=0; k<dim->ktot; k++){

    idx = j*dim->jstride + k*dim->kstride;
    
    ad_timestep_b(q[idx], qb[idx], grid->xy[idx], grid->xy[idx+jstride], grid->xy[idx+kstride],
    		  euler->inputs->cfl, &dt[idx], &dtb[idx]);

  }
  }

  for(j=0; j<dim->jtot; j++){
  for(k=0; k<dim->ktot; k++){
    idx = j*dim->jstride + k*dim->kstride;
    residual += qb[idx][0]*qb[idx][0];
    residual += qb[idx][1]*qb[idx][1];
    residual += qb[idx][2]*qb[idx][2];
    residual += qb[idx][3]*qb[idx][3];
    // add contribution from cost function
    qb[idx][0] = -qb2[idx][0] + qb[idx][0];
    qb[idx][1] = -qb2[idx][1] + qb[idx][1];
    qb[idx][2] = -qb2[idx][2] + qb[idx][2];
    qb[idx][3] = -qb2[idx][3] + qb[idx][3];
  }
  }

  this->boundary_conditions(false);

  // for(j=dim->nghost; j<dim->jmax+dim->nghost; j++){
  // for(k=dim->nghost; k<dim->kmax+dim->nghost; k++){
  for(j=0; j<dim->jtot; j++){
  for(k=0; k<dim->ktot; k++){

    idx = j*dim->jstride + k*dim->kstride;

    // update the residual
    rhsb[idx][0] = rhsb[idx][0] + qb[idx][0]*dt[idx];
    rhsb[idx][1] = rhsb[idx][1] + qb[idx][1]*dt[idx];
    rhsb[idx][2] = rhsb[idx][2] + qb[idx][2]*dt[idx];
    rhsb[idx][3] = rhsb[idx][3] + qb[idx][3]*dt[idx];
        
  }
  }
  
  this->boundary_conditions(false);
  
  step_number++;
  
  if(step_number % euler->inputs->resid == 0){
    residual = sqrt(residual/(dim->jtot*dim->ktot));
    printf("%4d : %15.6e\n", step_number, residual);
  }
  
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
