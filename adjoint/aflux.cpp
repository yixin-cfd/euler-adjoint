#include "adjoint.hpp"
#include "tmpfluxb.h"

// #define DO_SIMPLE
#define EPS 0.25

void adj_flux(double* q_l,   double* q_r, 
	      double* psi_l, double* psi_r, 
	      double* f_l,   double* f_r,
	      Dim *dim, double S[2], double A[4][4], int mini, int maxi){

  double Sx, Sy, phi, u, v, V, E;
  double a1;
  double a2 = GAMMA-1;
  double a3 = GAMMA-2;

  Sx         = S[0];
  Sy         = S[1];

  double dpsi0, dpsi1, dpsi2, dpsi3;
  dpsi0      = 0.5*(psi_l[0]*(1-mini) - psi_r[0]*(1-maxi));
  dpsi1      = 0.5*(psi_l[1]*(1-mini) - psi_r[1]*(1-maxi));
  dpsi2      = 0.5*(psi_l[2]*(1-mini) - psi_r[2]*(1-maxi));
  dpsi3      = 0.5*(psi_l[3]*(1-mini) - psi_r[3]*(1-maxi));

  // ---------------------------------------------------------------------
  // LEFT
  //
  u          = q_l[1]/q_l[0];
  v          = q_l[2]/q_l[0];
  E          = q_l[3]/q_l[0];

  V          = Sx*u + Sy*v;
  phi        = 0.5*a2*(u*u + v*v);
  a1         = GAMMA*E - phi;
  // build the Jacobian
  A[0][0]    = 0.0;
  A[0][1]    = Sx;
  A[0][2]    = Sy;
  A[0][3]    = 0.0;
  A[1][0]    = Sx*phi - u*V;
  A[1][1]    = V-a3*Sx*u;
  A[1][2]    = Sy*u-a2*Sx*v;
  A[1][3]    = a2*Sx;
  A[2][0]    = Sy*phi-v*V;
  A[2][1]    = Sx*v-a2*Sy*u;
  A[2][2]    = V-a3*Sy*v;
  A[2][3]    = a2*Sy;
  A[3][0]    = V*(phi-a1);
  A[3][1]    = Sx*a1-a2*u*V;
  A[3][2]    = Sy*a1-a2*v*V;
  A[3][3]    = GAMMA*V;

  // We want  [ A^T ]*[dpsi] - dissipation
  f_l[0]  = A[0][0]*dpsi0 + A[1][0]*dpsi1 + A[2][0]*dpsi2 + A[3][0]*dpsi3;
  f_l[1]  = A[0][1]*dpsi0 + A[1][1]*dpsi1 + A[2][1]*dpsi2 + A[3][1]*dpsi3;
  f_l[2]  = A[0][2]*dpsi0 + A[1][2]*dpsi1 + A[2][2]*dpsi2 + A[3][2]*dpsi3;
  f_l[3]  = A[0][3]*dpsi0 + A[1][3]*dpsi1 + A[2][3]*dpsi2 + A[3][3]*dpsi3;


  // ---------------------------------------------------------------------
  // RIGHT
  //
  u          = q_r[1]/q_r[0];
  v          = q_r[2]/q_r[0];
  E          = q_r[3]/q_r[0];

  V          = Sx*u + Sy*v;
  phi        = 0.5*a2*(u*u + v*v);
  a1         = GAMMA*E - phi;
  // build the Jacobian
  A[0][0]    = 0.0;
  A[0][1]    = Sx;
  A[0][2]    = Sy;
  A[0][3]    = 0.0;
  A[1][0]    = Sx*phi - u*V;
  A[1][1]    = V-a3*Sx*u;
  A[1][2]    = Sy*u-a2*Sx*v;
  A[1][3]    = a2*Sx;
  A[2][0]    = Sy*phi-v*V;
  A[2][1]    = Sx*v-a2*Sy*u;
  A[2][2]    = V-a3*Sy*v;
  A[2][3]    = a2*Sy;
  A[3][0]    = V*(phi-a1);
  A[3][1]    = Sx*a1-a2*u*V;
  A[3][2]    = Sy*a1-a2*v*V;
  A[3][3]    = GAMMA*V;

  f_r[0]  = A[0][0]*dpsi0 + A[1][0]*dpsi1 + A[2][0]*dpsi2 + A[3][0]*dpsi3;
  f_r[1]  = A[0][1]*dpsi0 + A[1][1]*dpsi1 + A[2][1]*dpsi2 + A[3][1]*dpsi3;
  f_r[2]  = A[0][2]*dpsi0 + A[1][2]*dpsi1 + A[2][2]*dpsi2 + A[3][2]*dpsi3;
  f_r[3]  = A[0][3]*dpsi0 + A[1][3]*dpsi1 + A[2][3]*dpsi2 + A[3][3]*dpsi3;

}

void adj_diss(double* q_l, double* q_r, double* psi_l, double* psi_r, double *d,
	      Dim *dim, double S[2]){

  double Sx, Sy, u, v, V, E;
  double a2 = GAMMA-1;
  double spec, c; // spectral radius, 

  Sx         = S[0];
  Sy         = S[1];
  u          = 0.5*(q_l[1]/q_l[0] + q_r[1]/q_r[0]);
  v          = 0.5*(q_l[2]/q_l[0] + q_r[2]/q_r[0]);
  E          = 0.5*(q_l[3]/q_l[0] + q_r[3]/q_r[0]);
  
  V          = Sx*u + Sy*v;
  c          = sqrt( a2*GAMMA*(E - 0.5*(u*u + v*v)) * (Sx*Sx + Sy*Sy) );
  spec       = std::abs(V)+c;

  d[0]       = 0.5*EPS*spec*(psi_r[0] - psi_l[0]);
  d[1]       = 0.5*EPS*spec*(psi_r[1] - psi_l[1]);
  d[2]       = 0.5*EPS*spec*(psi_r[2] - psi_l[2]);
  d[3]       = 0.5*EPS*spec*(psi_r[3] - psi_l[3]);

}


void Adjoint::aflux(){

  int j, k, idx, idx1, idx2;
  
  double f1[4], f2[4], d1[4], d2[4];

  double A[4][4];

  double fac1 = 0.5;
  double fac2 = 1.0;

  int mini, maxi;

  //
  // J-Direction
  //
  for(k=dim->nghost; k< dim->kmax+dim->nghost; k++){
  for(j=dim->nghost; j<=dim->jmax+dim->nghost; j++){

    idx  = j*dim->jstride + k*dim->kstride;
    idx1 = idx-dim->jstride;
    mini = (j == dim->nghost);
    maxi = (j == dim->nghost + dim->jmax);

    adj_flux(q[idx1], q[idx], psi[idx1], psi[idx], f1, f2, dim, grid->Sj[idx], A, mini, maxi);
    rhs[idx1][0] += f1[0];
    rhs[idx1][1] += f1[1];
    rhs[idx1][2] += f1[2];
    rhs[idx1][3] += f1[3];
    rhs[idx][0]  += f2[0];
    rhs[idx][1]  += f2[1];
    rhs[idx][2]  += f2[2];
    rhs[idx][3]  += f2[3];

  }
  }

  //
  // K-direction
  //
  for(j=dim->nghost; j< dim->jmax+dim->nghost; j++){
  for(k=dim->nghost; k<=dim->kmax+dim->nghost; k++){

    idx  = j*dim->jstride + k*dim->kstride;
    idx1 = idx-dim->kstride;
    mini = (k == dim->nghost);
    maxi = (k == dim->nghost + dim->kmax);

    adj_flux(q[idx1], q[idx], psi[idx1], psi[idx], f1, f2, dim, grid->Sk[idx], A, mini, maxi);
    rhs[idx1][0] += f1[0];
    rhs[idx1][1] += f1[1];
    rhs[idx1][2] += f1[2];
    rhs[idx1][3] += f1[3];
    rhs[idx][0]  += f2[0];
    rhs[idx][1]  += f2[1];
    rhs[idx][2]  += f2[2];
    rhs[idx][3]  += f2[3];

  }
  }

  //
  // Dissipation: J-Direction 
  //
  for(k=dim->nghost; k< dim->kmax+dim->nghost; k++){
  for(j=dim->nghost; j<=dim->jmax+dim->nghost; j++){

    idx  = j*dim->jstride + k*dim->kstride;
    idx1 = idx-dim->jstride;

    adj_diss(q[idx1], q[idx], psi[idx1], psi[idx], d1, dim, grid->Sj[idx]);
    for(int i=0; i<4; i++){
      rhs[idx ][i] += d1[i];
      rhs[idx1][i] -= d1[i];
    }
  }
  }
  //
  // Dissipation: K-direction
  //
  for(j=dim->nghost; j< dim->jmax+dim->nghost; j++){
  for(k=dim->nghost; k<=dim->kmax+dim->nghost; k++){

    idx  = j*dim->jstride + k*dim->kstride;
    idx1 = idx-dim->kstride;

    adj_diss(q[idx1], q[idx], psi[idx1], psi[idx], d1, dim, grid->Sk[idx]);
    for(int i=0; i<4; i++){
      rhs[idx ][i] += d1[i];
      rhs[idx1][i] -= d1[i];
    }
  }
  }

}
