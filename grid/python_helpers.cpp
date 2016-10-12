#include "python_helpers.hpp"
#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL flow_ARRAY_API
#include <numpy/ndarrayobject.h>


boost::python::object double_array_to_numpy(double* x, int len){

  npy_intp dims[1];
  dims[0] = len;

  PyObject *pyObj = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, (void*)x);

  boost::python::object bo(boost::python::handle<>(boost::python::borrowed(pyObj)));

  return bo;

}

boost::python::object int_array_to_numpy(int* x, int len){

  npy_intp dims[1];
  dims[0] = len;

  PyObject *pyObj = PyArray_SimpleNewFromData(1, dims, NPY_INT, (void*)x);

  boost::python::object bo(boost::python::handle<>(boost::python::borrowed(pyObj)));

  return bo;

}


void numpy_to_double_array(boost::python::object bo, double* x, int *len){

  PyArrayObject *arr;
  arr = (PyArrayObject*) bo.ptr();

  int ndim;
  npy_intp *dims;

  len[0] = 1;

  ndim = PyArray_NDIM(arr);
  dims = PyArray_DIMS(arr);

  for(int i=0; i<ndim; i++){
    len[0] *= dims[i];
  }

  x = (double*)PyArray_DATA(arr);

}
