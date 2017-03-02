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
    		   grid->xy[idx+dim->jstride], xyb[idx+dim->jstride],
    		   grid->xy[idx], xyb[idx], k, 0);
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

	// if(j == 5 && k == 1){
	//   printf("__ %20.16e %20.16e %20.16e %20.16e \n", 
	// 	 qb[idx][0], qb[idx][1], qb[idx][2], qb[idx][3]);
	// }


    }
    }

    //
    // K-direction
    //
    for(j=dim->nghost;   j< dim->jmax+dim->nghost; j++){
    for(k=dim->nghost;   k<=dim->kmax+dim->nghost; k++){

    	idx = j*dim->jstride + k*dim->kstride;

	// if(j == 1 && k == 1){
	// }


    	jkroeflux_b(q[idx-dim->kstride], qb[idx-dim->kstride], q[idx], qb[idx],
    		   rhs[idx-dim->kstride], rhsb[idx-dim->kstride],
    		   rhs[idx], rhsb[idx], dim,
    		   grid->xy[idx+dim->jstride], grid->xy[idx], k, 0);

	// if(j == 5 && k == 1){
	//   // printf("__ %20.16e %20.16e %20.16e %20.16e \n", 
	//   // 	 qb[idx][0], qb[idx][1], qb[idx][2], qb[idx][3]);
	//   printf("__ %20.16e %20.16e %20.16e %20.16e \n", 
	// 	 rhsb[idx][0], rhsb[idx][1], rhsb[idx][2], rhsb[idx][3]);
	// }



    }
    }
  }

}
