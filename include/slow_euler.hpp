#ifndef SLOW_EULER_H
#define SLOW_EULER_H
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <boost/python.hpp>
#include <stdio.h>
#include <string>
#include "grid.hpp"
#include "structures.h"

class Slow_Euler {

  void read_inputs(std::string yaml_string);
  void init();
  // void timestep();
  void boundary_conditions();
  void flux();
  double rhs_times_dt();
  void update_q();
  void dadi();

 public:
  Slow_Euler(Grid *g, std::string yaml_string);
  ~Slow_Euler();
  void say_hello();
  void write_solution(std::string name);
  void take_steps(int n);
  void go();
  void step();
  void save_restart(std::string fname);
  void read_restart(std::string fname);
  // int write_to_file(std::string s);
  // void metrics();

  boost::python::object pressure();
  boost::python::object Cl_Cd_Cm();

  int step_number;

  double (*f)[4];
  double (*q)[4];
  double (*rhs)[4];
  double *scratch, *dt;
  Dim *dim;
  BC *bc;
  Grid *grid;
  Inputs *inputs;

};

#endif
