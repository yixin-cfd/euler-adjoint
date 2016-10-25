#include "meshgen.hpp"
#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL mesh_ARRAY_API
#include <numpy/ndarrayobject.h>
#include "python_helpers.hpp"

int MeshGen::write_to_file(std::string s){

  FILE *fid;
  int pts = dim->pts;

  try {
    fid = fopen(s.c_str(),"w");
  } catch(int error){
    printf("*** Could not open file for writing\n");
    return 1;
  }

  int count = 0;
  fprintf(fid, "%d %d\n", dim->jtot, dim->ktot);

  for(int i=0; i<pts; i++){
    fprintf(fid, "%25.16e ", x[i]);
    count++;
    if(count%2 == 0)
      fprintf(fid, "\n");
  }
  for(int i=0; i<pts; i++){
    fprintf(fid, "%25.16e ", y[i]);
    count++;
    if(count%2 == 0)
      fprintf(fid, "\n");
  }

  fclose(fid);

  return 0;

}

boost::python::object MeshGen::get_mesh(){

  double (*xy)[2] = new double[dim->pts][2];

  int tmp_dims[] = {dim->ktot, dim->jtot, 2};

  for(int i=0; i<dim->pts; i++){
    xy[i][0] = x[i];
    xy[i][1] = y[i];
  }

  bool borrowed = false;
  boost::python::object bo;

  DOUBLE_TO_NUMPY(xy, bo, tmp_dims, 3, borrowed);

  return bo;

}

