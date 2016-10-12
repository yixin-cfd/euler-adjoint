#include "euler.hpp"

void Euler::take_steps(int n){
  for(int i=0; i<n; i++){
    step();
    step_number++;
  }
}

void Euler:: go(){
  this->take_steps(inputs->steps);
}

void Euler::step(){
  
  double residual;

  memset(rhs, 0, dim->pts*4*sizeof(double));

  this->timestep();

  this->boundary_conditions();

  this->flux();

  residual = this->rhs_times_dt();

  if(step_number % 1 == 0){

    if(residual != residual){
      throw 234;
    }

    printf("%d : %12.8e\n", step_number, residual);
  }

  this->dadi();

  this->update_q();
  

}
