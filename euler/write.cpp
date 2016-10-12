#include "euler.hpp"
#include <stdio.h>

#define DEBUG 1

void Euler::write_solution(std::string fname){

  int j, k, idx, var;
  FILE *fid;

  int jstride = dim->jstride;
  int kstride = dim->kstride;

  fid = fopen(fname.c_str(), "w");

  fprintf(fid,"VARIABLES = X,Y,RHO,RHO-U,RHO-V,E\n");
  fprintf(fid,"Zone i=%d,j=%d,F=BLOCK\n", dim->jtot, dim->ktot);
  fprintf(fid,"VARLOCATION = ([1-2]=NODAL,[3-6]=CELLCENTERED)\n");


  // X
  for (k=0;k<dim->ktot;k++){
    for (j=0;j<dim->jtot;j++){
      idx = j*jstride + k*kstride;
      fprintf(fid, "%25.16e\n", grid->xy[idx][0]);
    }
  }
  // Y
  for (k=0;k<dim->ktot;k++){
    for (j=0;j<dim->jtot;j++){
      idx = j*jstride + k*kstride;
      fprintf(fid, "%25.16e\n", grid->xy[idx][1]);
    }
  }
  int errc = 0;
  // all vars
  for (var=0;var<4;var++){
    for (k=0;k<dim->ktot-1;k++){
      for (j=0;j<dim->jtot-1;j++){
	idx = j*jstride + k*kstride;
	if(q[idx][var] != q[idx][var]){
	  // printf("found a nan at %d %d\n", j, k);
	  errc++;
	}
	if(errc < 5){
	  fprintf(fid, "%25.16e\n", q[idx][var]);	  
	} else {
	  printf("TOO MANY ERRORS, NOT FINISHING SOLUTION\n");
	  fclose(fid);
	  return;
	}
      }
    }
  }

  // for (k=0;k<dim->ktot-1;k++){
  //   for (j=0;j<dim->jtot-1;j++){
  //     idx = j*jstride + k*kstride;
  //     fprintf(fid, "%25.16e\n", scratch[idx]);
  //   }
  // }


  fclose(fid);

}
