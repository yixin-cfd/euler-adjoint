#ifndef ADADJ_H
#define ADADJ_H
#include <Python.h>
#include <boost/python.hpp>
#include <stdio.h>
#include <string>
#include "euler.hpp"
#include "structures.h"

class ADadj {

  Euler *euler;
  Dim *dim;
  Grid *grid;
  BC* wall;

  double (*f)[4];
  double (*q)[4];
  double (*rhs)[4];
  double *dt;

  double (*fb)[4];
  double (*qb)[4];
  double (*qb2)[4];
  double (*rhsb)[4];
  double (*xyb)[2];
  double *dtb;

  double *scratch;
  double *p_des;

  int step_number;

  void boundary_conditions(bool xbar);
  void flux(bool xbar);

public:
  ADadj(Euler *e);
  ~ADadj();
  void init(boost::python::object po);
  double  step();
  double check();
  void read_restart(std::string n);
  void save_restart(std::string n);
  void take_steps(int sn);

};


#endif
