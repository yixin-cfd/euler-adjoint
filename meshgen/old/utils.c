#include <stdio.h>
#include <math.h>
#include <sys/time.h>

double l2_norm_sq(double *vec, int len){
  int i;
  double norm = 0.0;
  for(i=0; i<len; i++){
    norm += vec[i]*vec[i];
  }
  return norm;
}

int stop_time(struct timeval *start, long *elapsed){
  struct timeval end;
  gettimeofday(&end, NULL);
  *elapsed = (end.tv_sec - (*start).tv_sec)*1000 + (end.tv_usec - (*start).tv_usec)/1000.0 + .5;
  return 0;
}

