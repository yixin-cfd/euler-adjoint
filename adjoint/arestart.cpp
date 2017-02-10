#include "adjoint.hpp"


void Adjoint::save_restart(std::string fname){

  int i, pts;
  FILE *fid;

  fid = fopen(fname.c_str(), "w");

  fprintf(fid, "%d %d\n", dim->jtot, dim->ktot);

  pts = dim->jtot*dim->ktot;

  for(i=0; i<pts;i++){

    fprintf(fid, "%24.16e %24.16e %24.16e %24.16e\n", 
	    psi[i][0], psi[i][1], psi[i][2], psi[i][3]);

  }

  fclose(fid);

}

void Adjoint::read_restart(std::string fname){
  int i, pts;
  FILE *fid;

  fid = fopen(fname.c_str(), "r");

  int test_jtot, test_ktot;
  fscanf(fid, "%d %d\n", &test_jtot, &test_ktot);

  if(test_jtot != dim->jtot || test_ktot != dim->ktot){
    printf("dimension mismatch, not reading restart\n");
    return;
  }

  pts = dim->jtot*dim->ktot;

  for(i=0; i<pts;i++){

    fscanf(fid, "%le %le %le %le\n", &psi[i][0], &psi[i][1], &psi[i][2], &psi[i][3]);

  }

  fclose(fid);

}
