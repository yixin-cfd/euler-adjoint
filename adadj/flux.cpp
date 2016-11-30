#include "adadj.hpp"
#include "adj_routines.h"

void ADadj::flux(bool xbar){
  
  int j, k, idx;

  if(xbar){
    //
    // J-direction
    //
    for(k=dim->nghost; k< dim->kmax+dim->nghost; k++){
      for(j=dim->nghost; j<=dim->jmax+dim->nghost; j++){

	idx = j*dim->jstride + k*dim->kstride;

	jroeflux_bx(q[idx-dim->jstride], qb[idx-dim->jstride], q[idx], qb[idx],
		    rhs[idx-dim->jstride], rhsb[idx-dim->jstride],
		    rhs[idx], rhsb[idx], dim,
		    grid->xy[idx], xyb[idx],
		    grid->xy[idx+dim->kstride], xyb[idx+dim->kstride], j);

      }
    }

    //
    // K-direction
    //
    for(j=dim->nghost; j< dim->jmax+dim->nghost; j++){
      for(k=dim->nghost; k<=dim->kmax+dim->nghost; k++){

	idx = j*dim->jstride + k*dim->kstride;

	kroeflux_bx(q[idx-dim->kstride], qb[idx-dim->kstride], q[idx], qb[idx],
		    rhs[idx-dim->kstride], rhsb[idx-dim->kstride],
		    rhs[idx], rhsb[idx], dim,
		    grid->xy[idx], xyb[idx],
		    grid->xy[idx+dim->jstride], xyb[idx+dim->jstride], k);

	// kroeflux(q[idx-dim->kstride], q[idx], rhs[idx-dim->kstride], rhs[idx], dim,
	// 	     grid->xy[idx], grid->xy[idx+dim->jstride], k);

      }
    }
  } else {
    //
    // J-direction
    //
    for(k=dim->nghost; k< dim->kmax+dim->nghost; k++){
      for(j=dim->nghost; j<=dim->jmax+dim->nghost; j++){

	idx = j*dim->jstride + k*dim->kstride;

	jroeflux_b(q[idx-dim->jstride], qb[idx-dim->jstride], q[idx], qb[idx],
		   rhs[idx-dim->jstride], rhsb[idx-dim->jstride],
		   rhs[idx], rhsb[idx], dim,
		   grid->xy[idx], grid->xy[idx+dim->kstride], j);

      }
    }

    //
    // K-direction
    //
    for(j=dim->nghost; j< dim->jmax+dim->nghost; j++){
      for(k=dim->nghost; k<=dim->kmax+dim->nghost; k++){

	idx = j*dim->jstride + k*dim->kstride;

	kroeflux_b(q[idx-dim->kstride], qb[idx-dim->kstride], q[idx], qb[idx],
		   rhs[idx-dim->kstride], rhsb[idx-dim->kstride],
		   rhs[idx], rhsb[idx], dim,
		   grid->xy[idx], grid->xy[idx+dim->jstride], k);

	// kroeflux(q[idx-dim->kstride], q[idx], rhs[idx-dim->kstride], rhs[idx], dim,
	// 	     grid->xy[idx], grid->xy[idx+dim->jstride], k);

      }
    }
  }

}
