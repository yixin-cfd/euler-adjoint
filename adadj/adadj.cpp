#include "adadj.hpp"
#include "euler.hpp"
#define PY_ARRAY_UNIQUE_SYMBOL adadj_ARRAY_API
#include <numpy/ndarrayobject.h>

ADadj::ADadj(Euler *e){

  import_array();

  this->euler = e;
  this->dim   = e->dim;
  this->grid  = e->grid;
  this->wall  = NULL;
  this->p_des = NULL;

  q       = new double[dim->pts][4];
  f       = new double[dim->pts][4];
  rhs     = new double[dim->pts][4];
  dt      = new double[dim->pts];

  qb      = new double[dim->pts][4];
  fb      = new double[dim->pts][4];
  rhsb    = new double[dim->pts][4];
  xyb     = new double[dim->pts][2];
  dtb     = new double[dim->pts];

  scratch = new double[dim->pts*4];

  printf("ADadj initialized\n");

}

ADadj::~ADadj(){

  delete f;
  delete q;
  delete rhs;
  delete dt;

  delete fb;
  delete qb;
  delete rhsb;
  delete dtb;

  delete scratch;

  if(p_des) delete p_des;

  printf("ADadj destroyed\n");
  
}

void ADadj::init(boost::python::object cp_desired_o){

  PyArrayObject *arr = (PyArrayObject*) cp_desired_o.ptr();

  int ndim;
  npy_intp *dims;

  int len = 1;

  ndim = PyArray_NDIM(arr);
  dims = PyArray_DIMS(arr);

  for(int i=0; i<ndim; i++){
    len *= dims[i];
  }

  double *cp_tmp = (double*)PyArray_DATA(arr);

  //
  // Now make sure we find the wall bc
  //
  bool found = false;
  for(int i=0; i<euler->inputs->nbc; i++){
    if(euler->bc[i].type == WALL_BC){
      found = true;
      wall = &euler->bc[i];
      break;
    }
  }
  if(not found){
    printf("No Wall Boundary Found\n");
    throw 431;
  }
  if(wall->face != KMIN_FACE){
    printf("Only KMIN faces can be pressure walls\n");
    throw 432;
  }

  int start, end;

  start = std::max(wall->js, dim->nghost);
  end   = std::min(wall->je, dim->jmax + dim->nghost - 1);

  int wall_pts = (end - start + 1);

  if(wall_pts != len){
    printf("You gave me the wrong number of pressure points!\n");
    throw 433;
  }

  p_des = new double[len];
  double dynp  = 0.5*euler->inputs->rho_inf*euler->inputs->M_inf*euler->inputs->M_inf;
  double p_inf = euler->inputs->p_inf;

  for(int i=0; i<len; i++){

    p_des[i] = cp_tmp[i]*dynp + p_inf;

  }
  
  printf("ok we got the pressure!\n");

}
