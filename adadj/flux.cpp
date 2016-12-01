#include "adadj.hpp"
#include "adj_routines.h"

void ADadj::flux(bool xbar){
  
  int j, k, idx;

  if(xbar){
    printf("this should not print\n");
    //
    // J-direction
    //
    for(k=dim->nghost; k< dim->kmax+dim->nghost; k++){
      for(j=dim->nghost; j<=dim->jmax+dim->nghost; j++){

	idx = j*dim->jstride + k*dim->kstride;

	jkroeflux_bx(q[idx-dim->jstride], qb[idx-dim->jstride], q[idx], qb[idx],
		     rhs[idx-dim->jstride], rhsb[idx-dim->jstride],
		     rhs[idx], rhsb[idx], dim,
		     grid->xy[idx], xyb[idx],
		     grid->xy[idx+dim->kstride], xyb[idx+dim->kstride], j, 1);

      }
    }

    //
    // K-direction
    //
    for(j=dim->nghost; j< dim->jmax+dim->nghost; j++){
      for(k=dim->nghost; k<=dim->kmax+dim->nghost; k++){

	idx = j*dim->jstride + k*dim->kstride;

	jkroeflux_bx(q[idx-dim->kstride], qb[idx-dim->kstride], q[idx], qb[idx],
		     rhs[idx-dim->kstride], rhsb[idx-dim->kstride],
		     rhs[idx], rhsb[idx], dim,
		     grid->xy[idx], xyb[idx],
		     grid->xy[idx+dim->jstride], xyb[idx+dim->jstride], k, 0);

      }
    }
  } else {
    //
    // J-direction
    //
    for(k=dim->nghost; k< dim->kmax+dim->nghost; k++){
    for(j=dim->nghost; j<=dim->jmax+dim->nghost; j++){

	idx = j*dim->jstride + k*dim->kstride;

	jkroeflux_b(q[idx-dim->jstride], qb[idx-dim->jstride], q[idx], qb[idx],
		    rhs[idx-dim->jstride], rhsb[idx-dim->jstride],
		    rhs[idx], rhsb[idx], dim,
		    grid->xy[idx], grid->xy[idx+dim->kstride], j, 1);

    }
    }

    // j = 8; k = 8;
    // idx = j*dim->jstride + k*dim->kstride;
    // printf("%15.8e %15.8e\n", qb[idx][0], qb[idx][1]);
    // printf("_ %15.8e %15.8e\n", qb2[idx][0], qb2[idx][1]);

    //
    // K-direction
    //
    for(j=dim->nghost;   j< dim->jmax+dim->nghost; j++){
    for(k=dim->nghost;   k<=dim->kmax+dim->nghost; k++){

    	idx = j*dim->jstride + k*dim->kstride;

    	jkroeflux_b(q[idx-dim->kstride], qb[idx-dim->kstride], q[idx], qb[idx],
    		   rhs[idx-dim->kstride], rhsb[idx-dim->kstride],
    		   rhs[idx], rhsb[idx], dim,
    		   grid->xy[idx+dim->jstride], grid->xy[idx], k, 0);

    }
    }
  }

}
