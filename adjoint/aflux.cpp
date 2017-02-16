#include "adjoint.hpp"

// #define DO_SIMPLE
#define EPS 0.125

void adj_flux(double* q_l, double* q_r, double* psi_l, double* psi_r, double* f, 
	      Dim *dim, double S[2], double A[4][4]){

  double Sx, Sy, phi, u, v, V, E;
  double a1;
  double a2 = GAMMA-1;
  double a3 = GAMMA-2;
  double dpsi0, dpsi1, dpsi2, dpsi3;
  double diss0, diss1, diss2, diss3;
  double spec, c; // spectral radius, 

  Sx         = S[0];
  Sy         = S[1];
  u          = 0.5*(q_l[1]/q_l[0] + q_r[1]/q_r[0]);
  v          = 0.5*(q_l[2]/q_l[0] + q_r[2]/q_r[0]);
  E          = 0.5*(q_l[3]/q_l[0] + q_r[3]/q_r[0]);
  
  V          = Sx*u + Sy*v;
  c          = sqrt( a2*GAMMA*(E - 0.5*(u*u + v*v)) * (Sx*Sx + Sy*Sy) );
  phi        = 0.5*a2*(u*u + v*v);
  a1         = GAMMA*E - phi;
  spec       = std::abs(V)+c;

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
	       
  dpsi0      = psi_r[0] - psi_l[0];
  dpsi1      = psi_r[1] - psi_l[1];
  dpsi2      = psi_r[2] - psi_l[2];
  dpsi3      = psi_r[3] - psi_l[3];

  diss0      = EPS*spec*dpsi0;
  diss1      = EPS*spec*dpsi1;
  diss2      = EPS*spec*dpsi2;
  diss3      = EPS*spec*dpsi3;

  // We want  [ A^T ]*[dpsi] - dissipation
  f[0]  = A[0][0]*dpsi0 + A[1][0]*dpsi1 + A[2][0]*dpsi2 + A[3][0]*dpsi3 - diss0;
  f[1]  = A[0][1]*dpsi0 + A[1][1]*dpsi1 + A[2][1]*dpsi2 + A[3][1]*dpsi3 - diss1;
  f[2]  = A[0][2]*dpsi0 + A[1][2]*dpsi1 + A[2][2]*dpsi2 + A[3][2]*dpsi3 - diss2;
  f[3]  = A[0][3]*dpsi0 + A[1][3]*dpsi1 + A[2][3]*dpsi2 + A[3][3]*dpsi3 - diss3;

}


void Adjoint::aflux(){

  int j, k, idx, idx1;

  double (*f)[4] = (double (*)[4])this->scratch;

  double A[4][4];

  //
  // J-Direction
  //
  for(k=dim->nghost; k< dim->kmax+dim->nghost; k++){
  for(j=dim->nghost; j<=dim->jmax+dim->nghost; j++){

    idx = j*dim->jstride + k*dim->kstride;
    idx1       = idx-dim->jstride;

    adj_flux(q[idx1], q[idx], psi[idx1], psi[idx], f[idx], dim, grid->Sj[idx], A);

    rhs[idx1][0] += (f[idx1][0] - f[idx][0])*(j>dim->nghost);
    rhs[idx1][1] += (f[idx1][1] - f[idx][1])*(j>dim->nghost);
    rhs[idx1][2] += (f[idx1][2] - f[idx][2])*(j>dim->nghost);
    rhs[idx1][3] += (f[idx1][3] - f[idx][3])*(j>dim->nghost);

  }
  }

  //
  // K-direction
  //
  for(j=dim->nghost; j< dim->jmax+dim->nghost; j++){
  for(k=dim->nghost; k<=dim->kmax+dim->nghost; k++){

    idx = j*dim->jstride + k*dim->kstride;
    idx1       = idx-dim->kstride;

    adj_flux(q[idx1], q[idx], psi[idx1], psi[idx], f[idx], dim, grid->Sk[idx], A);

    rhs[idx1][0] += (f[idx1][0] - f[idx][0])*(k>dim->nghost);
    rhs[idx1][1] += (f[idx1][1] - f[idx][1])*(k>dim->nghost);
    rhs[idx1][2] += (f[idx1][2] - f[idx][2])*(k>dim->nghost);
    rhs[idx1][3] += (f[idx1][3] - f[idx][3])*(k>dim->nghost);

  }
  }

}
