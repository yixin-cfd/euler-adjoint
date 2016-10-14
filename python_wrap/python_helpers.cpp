#include "python_helpers.hpp"
#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL flow_ARRAY_API
#include <numpy/ndarrayobject.h>

template <typename T>
boost::python::object array_to_numpy(T* x, int* dim, int ndim, bool borrowed){

  npy_intp dims[ndim];

  for(int i=0; i<ndim; i++){
    dims[i] = dim[i];
  }

  PyObject *pyObj;
  if(boost::is_same<T,int>::value){
    pyObj = PyArray_SimpleNewFromData(ndim, dims, NPY_INT, (void*)x);
  } else  if(boost::is_same<T,double>::value){
    pyObj = PyArray_SimpleNewFromData(ndim, dims, NPY_DOUBLE, (void*)x);
  } else {
    throw 234;
  }

  boost::python::object bo;

  if(borrowed){
    bo = boost::python::object(boost::python::handle<>(boost::python::borrowed(pyObj)));
  } else {
    PyArray_ENABLEFLAGS((PyArrayObject*)pyObj, NPY_ARRAY_OWNDATA);
    bo = boost::python::object(boost::python::handle<>(pyObj));
  } 

  return bo;
  
}

void compilerhelp(){
  array_to_numpy<double>(NULL, NULL, 1, true);
  array_to_numpy<int>(NULL, NULL, 1, true);
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
