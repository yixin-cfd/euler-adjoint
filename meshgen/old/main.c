#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include "headers.h"

#define ALGEBRAIC   4000
#define LAPLACE     5000
#define POISSON     6000

int main(int argc, char *argv[]){

  double *x, *y, *rhs;
  int jtot, ktot, jle, n, i, case_type;
  char garbage[LINE_MAX];
  char residfile[100];
  FILE *fid;
  double ds, xsf, res, w, xint;
  int nmax=400;

  fid = fopen("input_file.inp", "r");
  fgets(garbage, LINE_MAX, fid);                    // skip a line
  fscanf(fid, "%d %d\n", &jtot, &ktot);             // read these values
  fgets(garbage, LINE_MAX, fid);                    //
  fscanf(fid, "%lf %lf\n", &xsf, &ds);              // read these values
  fgets(garbage, LINE_MAX, fid);                    //
  fscanf(fid, "%d\n",&case_type);                   // case to run
  fgets(garbage, LINE_MAX, fid);                    //
  fscanf(fid, "%lf\n",&w);                          // omega for slor
  fgets(garbage, LINE_MAX, fid);                    //
  fscanf(fid, "%lf\n",&xint);                       // x_int
  fclose(fid);

  if(case_type == 0) case_type=ALGEBRAIC;
  if(case_type == 1) case_type=LAPLACE;
  if(case_type == 2) case_type=POISSON;

  printf("jtot, ktot: %d \t%d\n", jtot, ktot);             // read these values
  printf("xsf , ds  : %6.4f \t%6.4f\n", xsf, ds);              // read these values
  if(case_type == ALGEBRAIC)
      printf("case type : Algebraic\n");
  else if(case_type == LAPLACE)
      printf("case type : Laplace\n");
  else
      printf("case type : Poisson (Middlecoff-Thomas)\n");

  printf("omega     : %lf\n",w);                           // omega for slor
  printf("x_int     : %lf\n",xint);

  jle    = (jtot-1)/2; // for 0 indexing

  x = (double *)malloc( jtot*ktot*sizeof(double) );
  y = (double *)malloc( jtot*ktot*sizeof(double) );
  rhs = (double *)malloc( 2*jtot*ktot*sizeof(double) );

  if(case_type == ALGEBRAIC) nmax = 0; // don't do any iterations

  //
  // Clear the solution file
  //
  init_cgns(); 
  init_grid(x,y,jtot,ktot,jle,ds,xsf,xint);

  //write_grid(0,x,y,jtot,ktot);

  sprintf(residfile, "residual_w%1.2f.dat", w);

  fid = fopen(residfile, "w");

  // start timing
  struct timeval start;
  long elapsed;
  gettimeofday(&start, NULL);

  i = 0;
  for(n=1;n<=nmax;n++){
    //
    // Do Routine
    //
    slor(w,x,y,jtot,ktot,case_type==POISSON);
    //
    // Check Residual, calculate RHS
    //
    if((n-1)%(10)==0){
      residual(x,y,rhs,jtot,ktot,case_type==POISSON);
      res = l2_norm_sq(rhs, 2*jtot*ktot);
      res = sqrt(res/(2*jtot*ktot));
      printf("Residual: %16.10e\n", res);
      fprintf(fid, "%d %16.10e\n", n, res);
    }
  }
  // stop timing
  stop_time(&start, &elapsed);
  printf("\n *** Execution time: %ld milliseconds\n", elapsed);

  fclose(fid);

  write_grid(1,x,y,jtot,ktot);

  free(x);
  free(y);
  free(rhs);
  return 0;
}
