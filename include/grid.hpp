#ifndef GRID_H
#define GRID_H
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <boost/python.hpp>
#include <stdio.h>
#include <string>
#include "structures.h"

class Grid {

 public:
  Grid(boost::python::object o, int ghost);
  ~Grid();
  int write_to_file(std::string s);
  void metrics();

  Dim *dim;

  double (*xy)[2];
  double (*Sj)[2];
  double (*Sk)[2];
  double *V;

};

#endif
