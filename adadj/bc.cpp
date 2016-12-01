#include "adadj.hpp"
#include "adj_routines.h"

void rhsb_periodic_bc(double (*rhsb)[4], Dim *dim, BCface face, int j, int k){

  int idx, pidx;
  int pad = dim->nghost;

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

  rhsb[idx][0]   = rhsb[pidx][0];
  rhsb[idx][1]   = rhsb[pidx][1];
  rhsb[idx][2]   = rhsb[pidx][2];
  rhsb[idx][3]   = rhsb[pidx][3];

}


void rhsb_wall_bc(double (*rhsb)[4], Dim *dim, int j, int k){

  int idx = j*dim->jstride + k*dim->kstride;
  
  int up = dim->kstride;

  rhsb[idx][0] = 2*rhsb[idx+up][0] - rhsb[idx+up+up][0];
  rhsb[idx][1] = 2*rhsb[idx+up][1] - rhsb[idx+up+up][1];
  rhsb[idx][2] = 2*rhsb[idx+up][2] - rhsb[idx+up+up][2];
  rhsb[idx][3] = 2*rhsb[idx+up][3] - rhsb[idx+up+up][3];

}

void ADadj::boundary_conditions(bool xbar){

  BC bc;

  int i, idx, j, k;
  int wall_idx = -1;

  for(i=0; i<euler->inputs->nbc; i++){

    bc = euler->bc[i];

    for(k=bc.ks; k<=bc.ke; k++){
    for(j=bc.js; j<=bc.je; j++){
	
      idx = j*dim->jstride + k*dim->kstride;
      
      // periodic relies upon other values, pass full q pointer
      if(bc.type == PERIODIC_BC){
      	periodic_bc_b(q, qb, rhs, dim, bc.face, j, k);
      	rhsb_periodic_bc(rhsb, dim, bc.face, j, k);
      }
      
      // // dirichlet: set constant
      // if(bc.type == FARFIELD_BC) farfield_bc(q[idx], euler->inputs);	
      //

      if(bc.type == WALL_BC){
	// if(xbar)
	//   wall_bc_bx(q, qb, dim, grid->xy[idx], xyb[idx], 
	// 	     grid->xy[idx+dim->jstride], xyb[idx+dim->jstride], j, k);
	// else 
	wall_bc_b(q, qb, rhs, dim, grid->xy[idx], grid->xy[idx+dim->jstride], j, k);
	rhsb_wall_bc(rhsb, dim, j, k);
      }
      
    }
    }

  }


}
