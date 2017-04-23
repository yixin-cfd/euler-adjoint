#include "euler.hpp"
#include "yaml-cpp/yaml.h"
#define PY_ARRAY_UNIQUE_SYMBOL euler_ARRAY_API
#include <numpy/ndarrayobject.h>


Euler::Euler(Grid *g, std::string yaml_string){

  import_array(); // make sure numpy is loaded

  YAML::Node inputs = YAML::Load(yaml_string);

  grid        = g;

  this->dim    = new Dim[1];
  this->dim[0] = g->dim[0];

  step_number = 0;

  printf("Dimensions: (%d, %d) including %d ghosts\n", dim->jtot, dim->ktot, dim->nghost);

  this->inputs = NULL;
  this->bc     = NULL;

  f       = new double[dim->pts][4];
  q       = new double[dim->pts][4];
  rhs     = new double[dim->pts][4];
  scratch = new double[dim->pts*4];
  dt      = new double[dim->pts];

  memset(f,   0, dim->pts*4*sizeof(double));
  memset(q,   0, dim->pts*4*sizeof(double));
  memset(rhs, 0, dim->pts*4*sizeof(double));
  memset(dt,  0, dim->pts*1*sizeof(double));

  read_inputs(yaml_string);
  init();

  FILE *f = fopen("residual.dat", "w");
  fclose(f);
  
}


void Euler::say_hello(){
  printf("Euler says hello!\n");
}


Euler::~Euler(){

  // printf("___ Euler: cleaning up\n");
  
  if(inputs) delete inputs;
  if(bc)     delete bc;

  delete f;
  delete q;
  delete scratch;
  delete rhs;
  delete dt;

}
