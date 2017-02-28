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

void farfield_bc(BC bc, double (*psi)[4], double (*rhs)[4], Dim *dim){

  int j, k, idx;

  for(k=bc.ks; k<=bc.ke; k++){
  for(j=bc.js; j<=bc.je; j++){
  
    idx  = j*dim->jstride + k*dim->kstride;
    
    rhs[idx][0]  = 0.0; 
    rhs[idx][1]  = 0.0; 
    rhs[idx][2]  = 0.0; 
    rhs[idx][3]  = 0.0;

    psi[idx][0]  = 0.0; 
    psi[idx][1]  = 0.0; 
    psi[idx][2]  = 0.0; 
    psi[idx][3]  = 0.0;

  }
  }

}

void wall_bc(BC bc, double (*psi)[4], double (*q)[4], double (*rhs)[4], 
	     double (*S)[2], Dim *dim){

  int j, k, idx, midx, widx, up;
  double Sx, Sy, imag;
  double tmpb, tmpb0, tmpb1, tmpb2;
  int jstride = dim->jstride;
  int kstride = dim->kstride;

  up = dim->kstride;

  for(k=bc.ks; k<=bc.ke; k++){
  for(j=bc.js; j<=bc.je; j++){
  
    idx  = j*dim->jstride + k*dim->kstride;
    midx = j*jstride + (2*dim->nghost - k - 1)*kstride;
    widx = j*jstride + dim->nghost*kstride;
  
    Sx   = S[widx][0];
    Sy   = S[widx][1];

    imag = 1.0 / (Sx*Sx + Sy*Sy);

    Sx  *= imag;
    Sy  *= imag;

    tmpb = rhs[idx][3];
    rhs[idx][3] = 0.0;
    rhs[midx][3] = rhs[midx][3] + tmpb;
    tmpb0 = rhs[idx][2];
    rhs[idx][2] = 0.0;
    rhs[midx][2] = rhs[midx][2] + (1.0-Sy*Sy*2.0)*tmpb0;
    rhs[midx][1] = rhs[midx][1] - Sx*2.0*Sy*tmpb0;
    tmpb1 = rhs[idx][1];
    rhs[idx][1] = 0.0;
    rhs[midx][1] = rhs[midx][1] + (1.0-Sx*Sx*2.0)*tmpb1;
    rhs[midx][2] = rhs[midx][2] - Sx*2.0*Sy*tmpb1;
    tmpb2 = rhs[idx][0];
    rhs[idx][0] = 0.0;
    rhs[midx][0] = rhs[midx][0] + tmpb2;

    // extrapolate in k-direction beyond wall
    psi[idx][0] = 2*psi[idx+up][0] - psi[idx+up+up][0];
    psi[idx][1] = 2*psi[idx+up][1] - psi[idx+up+up][1];
    psi[idx][2] = 2*psi[idx+up][2] - psi[idx+up+up][2];
    psi[idx][3] = 2*psi[idx+up][3] - psi[idx+up+up][3];
  }
  }
    
}


void Adjoint::boundary_conditions(){

  for(int i=0; i<euler->inputs->nbc; i++){

    if (euler->bc[i].type == FARFIELD_BC){ 
      farfield_bc(euler->bc[i], psi, rhs, dim);
      continue;
    }

    switch(euler->bc[i].face){
    case JMIN_FACE:
      if (euler->bc[i].type == PERIODIC_BC) periodic_bc<JMIN_FACE>(euler->bc[i], psi, dim);
      break;
    case JMAX_FACE:
      if (euler->bc[i].type == PERIODIC_BC) periodic_bc<JMAX_FACE>(euler->bc[i], psi, dim);
      break;
    case KMIN_FACE:
      if (euler->bc[i].type == PERIODIC_BC) periodic_bc<KMIN_FACE>(euler->bc[i], psi, dim);
      if (euler->bc[i].type == WALL_BC)     wall_bc(euler->bc[i], psi, q, rhs, grid->Sk, dim);
      break;
    case KMAX_FACE:
      if (euler->bc[i].type == PERIODIC_BC) periodic_bc<KMAX_FACE>(euler->bc[i], psi, dim);
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
