#include "slow_euler.hpp"


void Slow_Euler::save_restart(std::string fname){

  int i, pts;
  FILE *fid;

  fid = fopen(fname.c_str(), "w");

  fprintf(fid, "%d %d\n", dim->jtot, dim->ktot);

  pts = dim->jtot*dim->ktot;

  for(i=0; i<pts;i++){

    fprintf(fid, "%24.16e %24.16e %24.16e %24.16e\n", 
	    q[i][0], q[i][1], q[i][2], q[i][3]);

  }

  fclose(fid);

}

void Slow_Euler::read_restart(std::string fname){
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

    fscanf(fid, "%le %le %le %le\n", &q[i][0], &q[i][1], &q[i][2], &q[i][3]);

  }

  fclose(fid);

}
