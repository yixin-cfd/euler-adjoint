#include "grid.hpp"
#define PY_ARRAY_UNIQUE_SYMBOL grid_ARRAY_API
#include <numpy/ndarrayobject.h>

// o is numpy array with shape [K, J, 2]
Grid::Grid(boost::python::object o, int ng){

  import_array();

  PyArrayObject *arr = (PyArrayObject*) o.ptr();

  int ndim, j, k;
  npy_intp *dims;

  ndim = PyArray_NDIM(arr);
  dims = PyArray_DIMS(arr);

  if(ndim != 3){
    printf("Numpy array must have 3 dimensions: x-y, j, k\n");
    throw 10;
  }
  if(dims[2] != 2){
    printf("Numpy array can only have 2 variables: x, y\n");
  }

  dim = new Dim[1];
  
  dim->nghost  = ng;
  dim->kmax    = dims[0]-1;
  dim->jmax    = dims[1]-1;
  dim->ktot    = dim->kmax + 2*dim->nghost;
  dim->jtot    = dim->jmax + 2*dim->nghost;

  dim->jstride = 1;
  dim->kstride = dim->jtot;

  int jstride = dim->jstride;
  int kstride = dim->kstride;

  dim->pts = dim->ktot * dim->jtot;

  double *in = (double*)PyArray_DATA(arr);

  xy = new double[dim->pts][2];
  Sj = new double[dim->pts][2];
  Sk = new double[dim->pts][2];
  V  = new double[dim->pts];

  int idx1, idx2;
  idx1 = 0;
  for(k=dim->nghost; k<=dim->ktot-dim->nghost; k++){
    for(j=dim->nghost; j<=dim->jtot-dim->nghost; j++){
      idx2 = j*jstride + k*kstride;
      xy[idx2][0] = in[idx1*2 + 0];
      xy[idx2][1] = in[idx1*2 + 1];
      
      idx1++;
    }
  }

  // handle j-ghost coords
  for(k=dim->nghost; k<=dim->ktot-dim->nghost; k++){
    for(j=dim->nghost-1; j>=0; j--){
      idx1 = dim->nghost*jstride + k*kstride;
      idx2 = j*jstride + k*kstride;
      xy[idx2][0] = 2*xy[idx2+jstride][0] - xy[idx2+2*jstride][0];
      xy[idx2][1] = 2*xy[idx2+jstride][1] - xy[idx2+2*jstride][1];
      // xy[idx2][0] = xy[idx1][0];
      // xy[idx2][1] = xy[idx1][1];
    }
    for(j=dim->jtot-dim->nghost+1; j<dim->jtot; j++){
      idx1 = (dim->jtot-dim->nghost)*jstride + k*kstride;
      idx2 = j*jstride + k*kstride;
      xy[idx2][0] = 2*xy[idx2-jstride][0] - xy[idx2-2*jstride][0];
      xy[idx2][1] = 2*xy[idx2-jstride][1] - xy[idx2-2*jstride][1];
      // xy[idx2][0] = xy[idx1][0];
      // xy[idx2][1] = xy[idx1][1];
    }    
  }


  // handle k-ghost coords
  for(j=0; j<dim->jtot; j++){
    for(k=dim->nghost-1; k>=0; k--){
      idx1 = j*jstride + dim->nghost*kstride;
      idx2 = j*jstride + k*kstride;
      xy[idx2][0] = 2*xy[idx2+kstride][0] - xy[idx2+2*kstride][0];
      xy[idx2][1] = 2*xy[idx2+kstride][1] - xy[idx2+2*kstride][1];
      // xy[idx2][0] = xy[idx1][0];
      // xy[idx2][1] = xy[idx1][1];
    }
    for(k=dim->ktot-dim->nghost+1; k<dim->ktot; k++){
      idx1 = j*jstride + (dim->ktot-dim->nghost)*kstride;
      idx2 = j*jstride + k*kstride;
      xy[idx2][0] = 2*xy[idx2-kstride][0] - xy[idx2-2*kstride][0];
      xy[idx2][1] = 2*xy[idx2-kstride][1] - xy[idx2-2*kstride][1];
      // xy[idx2][0] = xy[idx1][0];
      // xy[idx2][1] = xy[idx1][1];
    }
  }


  this->metrics();


}


Grid::~Grid(){
  
  // printf("___ Grid: deleting xy!!\n");
  delete[] xy;

}
