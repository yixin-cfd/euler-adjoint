#include "adjoint.hpp"

template<BCface face>
void periodic_bc(BC bc, double (*psi)[4], Dim *dim){

  int j, k, idx, pidx;
  int pad = dim->nghost;

  for(k=bc.ks; k<=bc.ke; k++){
  for(j=bc.js; j<=bc.je; j++){
      
      idx = j*dim->jstride + k*dim->kstride;

      // find the periodic idx
      if(face == JMIN_FACE){
      	pidx = (dim->jtot - 2*pad +j)*dim->jstride + k*dim->kstride;
      }
      if(face == JMAX_FACE){
      	pidx = (j - dim->jtot + 2*pad)*dim->jstride + k*dim->kstride;
      }
      if(face == KMIN_FACE){
      	pidx = j*dim->jstride + (dim->ktot - 2*pad + k)*dim->kstride;
      }
      if(face == KMAX_FACE){
      	pidx = j*dim->jstride + (k - dim->ktot + 2*pad)*dim->kstride;
      }

      psi[idx][0] = psi[pidx][0];
      psi[idx][1] = psi[pidx][1];
      psi[idx][2] = psi[pidx][2];
      psi[idx][3] = psi[pidx][3];

  }
  }
}

// template<BCface face>
// void periodic_bc(BC bc, double (*psi)[4], Dim *dim){

//   int j, k, idx, pidx;
//   int pad = dim->nghost;

//   for(k=bc.ks; k<=bc.ke; k++){
//   for(j=bc.js; j<=bc.je; j++){
      
//       idx = j*dim->jstride + k*dim->kstride;

//       // find the periodic idx
//       if(face == JMIN_FACE){
//       	pidx = (dim->jtot - 2*pad +j)*dim->jstride + k*dim->kstride;
//       }
//       if(face == JMAX_FACE){
//       	pidx = (j - dim->jtot + 2*pad)*dim->jstride + k*dim->kstride;
//       }
//       if(face == KMIN_FACE){
//       	pidx = j*dim->jstride + (dim->ktot - 2*pad + k)*dim->kstride;
//       }
//       if(face == KMAX_FACE){
//       	pidx = j*dim->jstride + (k - dim->ktot + 2*pad)*dim->kstride;
//       }

//       psi[idx][0] = psi[pidx][0];
//       psi[idx][1] = psi[pidx][1];
//       psi[idx][2] = psi[pidx][2];
//       psi[idx][3] = psi[pidx][3];

//   }
//   }
// }

void Adjoint::boundary_conditions(){

  for(int i=0; i<euler->inputs->nbc; i++){

    switch(euler->bc[i].face){
    case JMIN_FACE:
      if (euler->bc[i].type == PERIODIC_BC) periodic_bc<JMIN_FACE>(euler->bc[i], psi, dim);
      // if (euler->bc[i].type == FARFIELD_BC) farfield_bc(euler->bc[i], q, dim, inputs);
      break;
    case JMAX_FACE:
      if (euler->bc[i].type == PERIODIC_BC) periodic_bc<JMAX_FACE>(euler->bc[i], psi, dim);
      // if (euler->bc[i].type == FARFIELD_BC) farfield_bc(euler->bc[i], q, dim, inputs);
      break;
    case KMIN_FACE:
      if (euler->bc[i].type == PERIODIC_BC) periodic_bc<KMIN_FACE>(euler->bc[i], psi, dim);
      // if (euler->bc[i].type == FARFIELD_BC) farfield_bc(euler->bc[i], q, dim, inputs);
      // if (euler->bc[i].type == WALL_BC)     wall_bc<KMIN_FACE>(euler->bc[i], q, dim, grid);
      break;
    case KMAX_FACE:
      if (euler->bc[i].type == PERIODIC_BC) periodic_bc<KMAX_FACE>(euler->bc[i], psi, dim);
      // if (euler->bc[i].type == FARFIELD_BC) farfield_bc(euler->bc[i], q, dim, inputs);
      // if (euler->bc[i].type == WALL_BC)     wall_bc<KMAX_FACE>(euler->bc[i], q, dim, grid);
      break;

    }

    // // slow start
    // if(ramp && euler->bc[i].type == WALL_BC){
    //   ratio = step_number * 1.0 / 30.0;
    //   ratio = (10.0 - 15.0*ratio + 6.0*ratio*ratio)*ratio*ratio*ratio;
    //   ramp_bc(ratio, euler->bc[i], q, dim, inputs);
    // }
    
  }


}
