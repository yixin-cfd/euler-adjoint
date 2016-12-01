#include "adadj.hpp"

void ADadj::save_restart(std::string fname){

  int i, pts;
  FILE *fid;

  fid = fopen(fname.c_str(), "w");

  fprintf(fid, "%d %d\n", dim->jtot, dim->ktot);

  pts = dim->jtot*dim->ktot;

  for(i=0; i<pts;i++){

    fprintf(fid, "%24.16e %24.16e %24.16e %24.16e\n", 
	    rhsb[i][0], rhsb[i][1], rhsb[i][2], rhsb[i][3]);

  }

  fclose(fid);

}

void ADadj::read_restart(std::string fname){
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

    fscanf(fid, "%le %le %le %le\n", &rhsb[i][0], &rhsb[i][1], &rhsb[i][2], &rhsb[i][3]);

  }

  fclose(fid);

}
