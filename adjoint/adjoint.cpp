#include "euler.hpp"
#include "yaml-cpp/yaml.h"
#define PY_ARRAY_UNIQUE_SYMBOL euler_ARRAY_API
#include <numpy/ndarrayobject.h>


Adjoint::Adjoint(Euler *e){

  import_array(); // make sure numpy is loaded

  this->euler = e;
  this->dim   = e->dim;
  this->grid  = e->grid;
  this->wall  = NULL;
  this->p_des = NULL;

  q       = new double[dim->pts][4];

  dt      = new double[dim->pts];
  psi     = new double[dim->pts][4];
  rhs     = new double[dim->pts][4];

  scratch = new double[dim->pts*4];

  step_number = 0;

  // use the Q from the Euler class
  memcpy(q, euler->q, dim->pts*4*sizeof(double));

  // zero everything else
  memset(rhsb, 0, 4*dim->pts*sizeof(double));


  printf("ADadj initialized\n");
  
}


void Adjoint::say_hello(){
  printf("Adjoint says hello!\n");
}


Adjoint::~Adjoint(){

  if(inputs) delete inputs;
  if(bc)     delete bc;

  delete f;
  delete q;
  delete scratch;
  delete rhs;
  delete dt;

}
