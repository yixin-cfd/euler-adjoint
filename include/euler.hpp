#ifndef EULER_H
#define EULER_H
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <boost/python.hpp>
#include <stdio.h>
#include <string>
#include "grid.hpp"
#include "structures.h"

class Euler {

  void read_inputs(std::string yaml_string);
  void init();
  void timestep();
  void boundary_conditions();
  void flux();
  double rhs_times_dt();
  void update_q();
  void dadi();

  Inputs *inputs;
  BC *bc;
  Grid *grid;
  Dim *dim;

 public:
  Euler(Grid *g, std::string yaml_string);
  ~Euler();
  void say_hello();
  void write_solution(std::string name);
  void take_steps(int n);
  void go();
  void step();
  // int write_to_file(std::string s);
  // void metrics();

  int step_number;

  double (*f)[4];
  double (*q)[4];
  double (*rhs)[4];
  double *scratch, *dt;

};

#endif
