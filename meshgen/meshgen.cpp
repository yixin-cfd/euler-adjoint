#include "meshgen.hpp"
#define PY_ARRAY_UNIQUE_SYMBOL mesh_ARRAY_API
#include <numpy/ndarrayobject.h>
// #include "python_helpers.hpp"

MeshGen::MeshGen(boost::python::object o, int ktot, double dist){

  import_array();

  PyArrayObject *arr = (PyArrayObject*) o.ptr();

  int ndim, j, k;
  npy_intp *dims;

  ndim = PyArray_NDIM(arr);
  dims = PyArray_DIMS(arr);

  if(ndim != 2){
    printf("Numpy array must have 2 dimensions: x-y, j\n");
    throw 10;
  }
  if(dims[1] != 2){
    printf("Numpy array can only have 2 variables: x, y\n");
  }

  dim = new Dim[1];
  
  dim->nghost  = 0;
  dim->jmax    = dims[0];
  dim->kmax    = ktot;
  dim->ktot    = dim->kmax;
  dim->jtot    = dim->jmax;
  dim->pts     = dim->ktot * dim->jtot;
  dim->jstride = 1;
  dim->kstride = dim->jtot;

  int jstride = dim->jstride;
  int kstride = dim->kstride;

  double *in = (double*)PyArray_DATA(arr);

  x   = new double[dim->pts];
  y   = new double[dim->pts];
  a   = new double[dim->jtot];
  b   = new double[dim->jtot];
  c   = new double[dim->jtot];
  d   = new double[dim->jtot];
  P   = new double[dim->pts];
  Q   = new double[dim->pts];
  rhs = new double[dim->pts][2];

  int idx1;
  k = 0;
  for(j=0; j<=dim->jtot; j++){
    idx1 = j*jstride + k*kstride;
    x[idx1] = in[j*2 + 0];
    y[idx1] = in[j*2 + 1];
  }

  this->init(dist);

}

MeshGen::~MeshGen(){

  delete x;
  delete y;

  delete a;
  delete b;
  delete c;
  delete d;
  delete P;
  delete Q;

  delete rhs;

}
