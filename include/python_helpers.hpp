#ifndef GRID_PYTHON_HELPER_H
#define GRID_PYTHON_HELPER_H
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <stdio.h>
#include <boost/python.hpp>
#include "boost/python/extract.hpp"
// http://docs.scipy.org/doc/numpy/reference/c-api.array.html#miscellaneous


template <typename T>
boost::python::object array_to_numpy(T* x, int* dim, int ndim, bool borrowed);

// boost::python::object double_array_to_numpy(double* x, int len);

boost::python::object int_array_to_numpy(int* x, int len);

void numpy_to_double_array(boost::python::object bo, double* x, int *len);
#endif
