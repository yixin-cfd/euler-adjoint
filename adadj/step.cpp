#include "adadj.hpp"
#include "adj_routines.h"
#include "euler_routines.h"
#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL adadj_ARRAY_API
#include <numpy/ndarrayobject.h>

void ADadj::take_steps(int nsteps){

  FILE *file;
  int i, j, k, idx;

  memset( qb2, 0, 4*dim->pts*sizeof(double));

  //
  // Derivative of cost function J wrt J is 1!
  //
  double dummy, liftd, residual, Jb = 1.0;

  double upx, upy;

  // figure out which direction is "up" for lift
  upy =  cos(euler->inputs->aoa*M_PI/180.0);
  upx = -sin(euler->inputs->aoa*M_PI/180.0);

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

    pressure_cost_b( q[idx], qb2[idx], dim, &dummy, &Jb, p_des[pj], euler->inputs);
    // lift_cost_b(q[idx], qb2[idx], grid->xy[idx], grid->xy[idx+dim->jstride], dim, &dummy, &Jb,
    // 		upx, upy);
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

    if(step_number % euler->inputs->resid == 0 || step_number == nsteps){
      residual = sqrt(residual/(dim->jtot*dim->ktot));
      liftd    = this->check();
      printf("%5d  %15.6e  %15.8e\n", step_number, residual, liftd);
      file = fopen("res_adj.dat", "a");
      fprintf(file, "%5d  %15.6e  %15.8e\n", step_number, residual, liftd);
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

  // for(j=0; j<dim->jtot-1; j++){
  // for(k=0; k<dim->ktot-1; k++){

  //   idx = j*dim->jstride + k*dim->kstride;
    
  //   ad_timestep_b(q[idx], qb[idx], grid->xy[idx], grid->xy[idx+jstride], grid->xy[idx+kstride],
  //   		  euler->inputs->cfl, &dt[idx], &dtb[idx]);

  // }
  // }

  this->boundary_conditions(false);

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
    qb[idx][0] = qb2[idx][0] + qb[idx][0];
    qb[idx][1] = qb2[idx][1] + qb[idx][1];
    qb[idx][2] = qb2[idx][2] + qb[idx][2];
    qb[idx][3] = qb2[idx][3] + qb[idx][3];

    // update the residual
    rhsb[idx][0] = rhsb[idx][0] - qb[idx][0]*dt[idx];
    rhsb[idx][1] = rhsb[idx][1] - qb[idx][1]*dt[idx];
    rhsb[idx][2] = rhsb[idx][2] - qb[idx][2]*dt[idx];
    rhsb[idx][3] = rhsb[idx][3] - qb[idx][3]*dt[idx];

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

double ADadj::check(){

  int i, j, k, idx;
  int jstride = dim->jstride;
  int kstride = dim->kstride;
  
  memset( qb2, 0, 4*dim->pts*sizeof(double));
  memset(  qb, 0, 4*dim->pts*sizeof(double));
  memset( dtb, 0,   dim->pts*sizeof(double));
  memset( xyb, 0, 2*dim->pts*sizeof(double));

  double (*xd)[2] = new double[dim->pts][2];
  double dalpha   = atan(1.0)/45.0;

  for(i=0; i<dim->pts; i++){
    xd[i][0] =  grid->xy[i][1];//*dalpha;
    xd[i][1] = -grid->xy[i][0];//*dalpha;
  }

  //
  // Derivative of cost function J wrt J is 1!
  //
  double dummy, residual, Jb = 1.0; double actual_cost=0.0;
  double upx, upy;

  // figure out which direction is "up" for lift
  upy =  cos(euler->inputs->aoa*M_PI/180.0);
  upx = -sin(euler->inputs->aoa*M_PI/180.0);

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

    pressure_cost_b( q[idx], qb2[idx], dim, &dummy, &Jb, p_des[pj], euler->inputs);
    // lift_cost_bx(q[idx], qb2[idx], grid->xy[idx], xyb[idx], 
    // 		 grid->xy[idx+dim->jstride], xyb[idx+dim->jstride], dim, &dummy, &Jb, upx, upy);
    // lift_cost(q[idx], grid->xy[idx], grid->xy[idx+dim->jstride], dim, &actual_cost, upx, upy);

  }

  // printf("lift is %lf\n", actual_cost*8);
  //printf("testing123 %e", actual_cost);

  this->flux(true);
  this->boundary_conditions(true);

  double liftd = 0.0;
  for(i=0; i<dim->pts; i++){

    liftd += xyb[i][0]*xd[i][0] + xyb[i][1]*xd[i][1];

  }

  delete xd;

  return liftd;
  
}

double ADadj::sens_xd(boost::python::object xdo){

  PyArrayObject *arr = (PyArrayObject*) xdo.ptr();

  int ndim, len;
  npy_intp *dims;

  ndim = PyArray_NDIM(arr);
  dims = PyArray_DIMS(arr);

  if(ndim != 3){
    printf("ndim not 3\n");
  }

  if(dims[2] != 2){
    printf("last dim should be 2!\n");
  }

  int xjstride = 1;
  int xkstride = dims[1];

  double (*xd)[2] = (double (*)[2])PyArray_DATA(arr);

  // ------------------------------------------------------------
  //
  // Now tha actual adjoint calculation:
  //
  int i, j, k, idx, xidx;
  int jstride = dim->jstride;
  int kstride = dim->kstride;
  
  memset( qb2, 0, 4*dim->pts*sizeof(double));
  memset(  qb, 0, 4*dim->pts*sizeof(double));
  memset( dtb, 0,   dim->pts*sizeof(double));
  memset( xyb, 0, 2*dim->pts*sizeof(double));

  //
  // Derivative of cost function J wrt J is 1!
  //
  double dummy, Jb = 1.0;
  double upx, upy, actual_cost;

  // figure out which direction is "up" for lift
  upy =  cos(euler->inputs->aoa*M_PI/180.0);
  upx = -sin(euler->inputs->aoa*M_PI/180.0);

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

    // pressure_cost(q[idx], dim, &actual_cost, p_des[pj], euler->inputs);

    pressure_cost_b( q[idx], qb2[idx], dim, &dummy, &Jb, p_des[pj], euler->inputs);
    // lift_cost_bx(q[idx], qb2[idx], grid->xy[idx], xyb[idx], 
    // 		 grid->xy[idx+dim->jstride], xyb[idx+dim->jstride], dim, &dummy, &Jb, upx, upy);

  }

  //printf("initial actual cost is %lf\n", actual_cost);
  //printf("testing123 %e", actual_cost);

  this->flux(true);
  this->boundary_conditions(true);

  double costd = 0.0;

  for(k=dim->nghost; k<=dim->kmax+dim->nghost; k++){
  for(j=dim->nghost; j<=dim->jmax+dim->nghost; j++){
    idx  = j*jstride + k*kstride;
    xidx = (j-dim->nghost)*xjstride + (k-dim->nghost)*xkstride;

    costd += xyb[idx][0]*xd[xidx][0] + xyb[idx][1]*xd[xidx][1];

  }
  }

  return costd;

}
