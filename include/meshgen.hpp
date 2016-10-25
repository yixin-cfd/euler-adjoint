#ifndef MESHGEN_H
#define MESHGEN_H
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <boost/python.hpp>
#include <stdio.h>
#include <string>
#include "structures.h"

class MeshGen {

  double *x, *y;
  double *a, *b, *c, *d, *P, *Q;
  double (*rhs)[2];
  Dim *dim;
  
 public:
  MeshGen(boost::python::object o, int ktot, double dist);
  ~MeshGen();
  void init(double dist);
  int write_to_file(std::string s);
  boost::python::object get_mesh();
  void poisson(int n);
};


#endif
