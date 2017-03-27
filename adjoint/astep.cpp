#include "adjoint.hpp"
#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL adjoint_ARRAY_API
#include <numpy/ndarrayobject.h>

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

    if(step_number % euler->inputs->resid == 0 || step_number == nsteps){
      residual = sqrt(residual/(dim->jtot*dim->ktot));
      liftd    = this->check();
      printf("%5d  %15.6e  %15.8e\n", step_number, residual, liftd);
      file = fopen("res_adj.dat", "a");
      fprintf(file, "%5d  %15.6e  %15.8e\n", step_number, residual, liftd);
      fclose(file);
      if(residual > 10) return;
    }

    // printf("Residual is %25.16e\n", residual);
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

  this->boundary_conditions();

  for(j=0; j<dim->jtot; j++){
  for(k=0; k<dim->ktot; k++){

    idx = j*dim->jstride + k*dim->kstride;

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
}


double Adjoint::sens_xd(boost::python::object xdo){

  PyArrayObject *arr = (PyArrayObject*) xdo.ptr();

  int ndim, j, k, idx, xidx;
  npy_intp *dims;
  double costd;

  ndim = PyArray_NDIM(arr);
  dims = PyArray_DIMS(arr);

  int xjstride = 1;
  int xkstride = dims[1];

  if(ndim != 3){
    printf("ndim not 3\n");
  }

  if(dims[2] != 2){
    printf("last dim should be 2!\n");
  }

  double (*xd)[2] = (double (*)[2])PyArray_DATA(arr);

  // get the deriv w.r.t. X and Y
  this->check();

  costd = 0.0;
  for(k=dim->nghost; k<=dim->kmax+dim->nghost; k++){
    for(j=dim->nghost; j<=dim->jmax+dim->nghost; j++){
      
      idx  = j*dim->jstride + k*dim->kstride;
      
      // the provided delta does NOT have ghosts
      xidx = (j-dim->nghost)*xjstride + (k-dim->nghost)*xkstride;
      
      costd += xyb[idx][0]*xd[xidx][0] + xyb[idx][1]*xd[xidx][1];

      if(j == 5 && k == 1){
      	printf("%d %d: %20.14e %20.14e \n", j, k, xyb[idx][0], xyb[idx][1]);
      }
      
    }
  }

  return costd;

}
