#include "structures.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPS 0.25

void tmpflux(double* q_l, double* q_r, double *rhs_m1, double* rhs, 
	     Dim *dim, double S[2], int mini, int maxi){
  double f[4];
  
// void tmpflux(double* q_l, double* q_r, double *f, Dim *dim, double S[2]){

  double rho_l  = q_l[0];
  double rho_r  = q_r[0];
  double inv_rho_l  = 1.0/rho_l;
  double inv_rho_r  = 1.0/rho_r;

  double u_l  = q_l[1]*inv_rho_l;
  double u_r  = q_r[1]*inv_rho_r;

  double v_l  = q_l[2]*inv_rho_l;
  double v_r  = q_r[2]*inv_rho_r;

  double e_l  = q_l[3];
  double e_r  = q_r[3];

  double p_l  = (GAMMA-1.0)*(e_l -  0.5*rho_l*(u_l*u_l + v_l*v_l));
  double p_r  = (GAMMA-1.0)*(e_r -  0.5*rho_r*(u_r*u_r + v_r*v_r));

  // face size and normals
  double mag = sqrt(S[0]*S[0] + S[1]*S[1]);
  double imag = 1.0/mag;
  double r1 = S[0]*imag;
  double r2 = S[1]*imag;

  double V_l = u_l *r1 + v_l *r2;
  double V_r = u_r *r1 + v_r *r2;

  double dF0=0.0, dF1=0.0, dF2=0.0, dF3=0.0;

  double c_l  = sqrt((GAMMA*p_l)*inv_rho_l);
  double c_r  = sqrt((GAMMA*p_r)*inv_rho_r);
  double eig_l, eig_r;
  eig_l =  sqrt(V_l*V_l) + c_l;
  eig_r =  sqrt(V_r*V_r) + c_r;

  double rad = 0.5*(eig_l + eig_r);

  dF0 = EPS * rad * (q_r[0] - q_l[0]);
  dF1 = EPS * rad * (q_r[1] - q_l[1]);
  dF2 = EPS * rad * (q_r[2] - q_l[2]);
  dF3 = EPS * rad * (q_r[3] - q_l[3]);

  double half_face = 0.5*mag;

  f[0] = half_face*( (rho_l*V_l             ) + (rho_r*V_r             ) - dF0);
  f[1] = half_face*( (rho_l*u_l*V_l + p_l*r1) + (rho_r*u_r*V_r + p_r*r1) - dF1);
  f[2] = half_face*( (rho_l*v_l*V_l + p_l*r2) + (rho_r*v_r*V_r + p_r*r2) - dF2);
  f[3] = half_face*( (e_l+p_l)*V_l + (e_r+p_r)*V_r                       - dF3);

  rhs_m1[0] = rhs_m1[0] + f[0]*(1-mini);
  rhs_m1[1] = rhs_m1[1] + f[1]*(1-mini);
  rhs_m1[2] = rhs_m1[2] + f[2]*(1-mini);
  rhs_m1[3] = rhs_m1[3] + f[3]*(1-mini);

  rhs[0]    = rhs[0]    - f[0]*(1-maxi);
  rhs[1]    = rhs[1]    - f[1]*(1-maxi);
  rhs[2]    = rhs[2]    - f[2]*(1-maxi);
  rhs[3]    = rhs[3]    - f[3]*(1-maxi);

  // if(j_or_k > dim->nghost){
  // }

  // if(is_j==1 && j_or_k < dim->nghost + dim->jmax){
  //   rhs[0]    = rhs[0] - f[0];
  //   rhs[1]    = rhs[1] - f[1];
  //   rhs[2]    = rhs[2] - f[2];
  //   rhs[3]    = rhs[3] - f[3];
  // }

  // if(is_j==0 && j_or_k < dim->nghost + dim->kmax){
  //   rhs[0]    = rhs[0] - f[0];
  //   rhs[1]    = rhs[1] - f[1];
  //   rhs[2]    = rhs[2] - f[2];
  //   rhs[3]    = rhs[3] - f[3];
  // }

}
