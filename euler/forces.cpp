#include "python_helpers.hpp"
#include "euler.hpp"

boost::python::object Euler::pressure(){

  BC wallbc;
  bool found = false;
  for(int i=0; i<inputs->nbc; i++){
    if(bc[i].type == WALL_BC){
      found = true;
      wallbc = bc[i];
      break;
    }
  }
  if(not found){
    return boost::python::object(); // None object
  }
  if(not wallbc.face == KMIN_FACE){
    printf("Only KMIN faces can be pressure walls\n");
    throw 432;
  }

  int start, end;

  start = std::max(wallbc.js, dim->nghost);
  end   = std::min(wallbc.je, dim->jmax + dim->nghost - 1);

  int wall_pts = (end - start + 1);

  // store x-coordinate and pressure
  double (*p)[2] = new double[wall_pts][2];

  int j, k, idx, idx1;
  k = dim->nghost;

  double p1, p2;
  double dynp = 0.5*inputs->rho_inf*inputs->M_inf*inputs->M_inf;

  for(j=start; j<=end; j++){
    
    idx = j*dim->jstride + k*dim->kstride;

    // x-coordinate
    idx1 = idx + dim->jstride;
    p[j-start][0] = 0.5*(grid->xy[idx][0] + grid->xy[idx1][0]);

    // pressure
    idx1 = idx - dim->kstride;
    p1 = (GAMMA-1.0)*(q[idx ][3] - 0.5*(q[idx ][1]*q[idx][1] + q[idx ][2]*q[idx ][2])/q[idx ][0]);
    p2 = (GAMMA-1.0)*(q[idx1][3] - 0.5*(q[idx1][1]*q[idx][1] + q[idx1][2]*q[idx1][2])/q[idx1][0]);
    p[j-start][1] = (0.5*(p1 + p2) - inputs->p_inf)/dynp;

  }

  // double (*x)[2] = new double[16][2];

  // for(int i=0; i<16; i++){
  //   for(int ii=0; ii<2; ii++){
  //     x[i][ii] = 1.0*i*ii;
  //   }
  // }

  int dims[] = {wall_pts, 2};

  // we dont want python to borrow the array, we want it to own the array
  bool borrowed = false; 
  boost::python::object bo = array_to_numpy<double>((double*)p, dims, 2, borrowed);

  return bo;

}

