#include "adadj.hpp"
#include "adj_routines.h"

void ADadj::boundary_conditions(bool xbar){

  BC bc;

  int i, idx, j, k;
  
  for(i=0; i<euler->inputs->nbc; i++){

    bc = euler->bc[i];

    for(k=bc.ks; k<=bc.ke; k++){
      for(j=bc.js; j<=bc.je; j++){
	
	idx = j*dim->jstride + k*dim->kstride;
	
	// periodic relies upon other values, pass full q pointer
	if(bc.type == PERIODIC_BC) periodic_bc_b(q, qb, dim, JMIN_FACE, j, k);
	
	// // dirichlet: set constant
	// if(bc.type == FARFIELD_BC) farfield_bc(q[idx], euler->inputs);
	
	// 
	if(bc.type == WALL_BC and xbar){
	  wall_bc_bx(q, qb, dim, grid->xy[idx], xyb[idx], 
		     grid->xy[idx+dim->jstride], xyb[idx+dim->jstride], j, k);
	} else if(bc.type == WALL_BC){
	  wall_bc_b(q, qb, dim, grid->xy[idx], grid->xy[idx+dim->jstride], j, k);
	}
      }
    }
    
  }

}
