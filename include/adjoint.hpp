#ifndef ADJOINT_H
#define ADJOINT_H
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <boost/python.hpp>
#include <stdio.h>
#include <string>
#include "euler.hpp"
#include "structures.h"

class Adjoint {

  Euler *euler;
  Dim *dim;
  Grid *grid;
  BC* wall;

  double *scratch, *dt;
  double *p_des;
  
  // void timestep();
  void boundary_conditions();
  void aflux();
  void adiss();
  double arhs_times_dt();
  void update_psi();
  void timestep();
  double check();
  void cost(double (*r)[4]);

 public:
  Adjoint(Euler *e);
  ~Adjoint();
  void init(boost::python::object po);
  void say_hello();
  void write_solution(std::string name);
  void take_steps(int n);
  void go();
  double step();
  void save_restart(std::string fname);
  void read_restart(std::string fname);
  // int write_to_file(std::string s);
  // void metrics();

  // boost::python::object pressure();
  // boost::python::object Cl_Cd_Cm();

  int step_number;

  double (*psi)[4];
  double (*q)[4];
  double (*rhs)[4];
  double (*rhs0)[4];
  double (*xyb)[2];

};

#endif
