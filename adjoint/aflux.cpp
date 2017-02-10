#include "adjoint.hpp"

// #define DO_SIMPLE
#define EPS 0.125

void flux(double* q_l, double* q_r, double* f, Dim *dim, double S[2]){

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

  double c_l  = sqrt((GAMMA*p_l)*inv_rho_l);
  double c_r  = sqrt((GAMMA*p_r)*inv_rho_r);

  // face size and normals
  double mag = sqrt(S[0]*S[0] + S[1]*S[1]);
  double imag = 1.0/mag;
  double r1 = S[0]*imag;
  double r2 = S[1]*imag;

  double V_l = u_l *r1 + v_l *r2;
  double V_r = u_r *r1 + v_r *r2;

  double eig_l  = std::abs(V_l) + c_l;
  double eig_r  = std::abs(V_r) + c_r;

  double rad = 0.5*(eig_l + eig_r);
  double dF0, dF1, dF2, dF3; // delta F, blazek eqn. 4.89-4.91

  dF0 = EPS * rad * (q_r[0] - q_l[0]);
  dF1 = EPS * rad * (q_r[1] - q_l[1]);
  dF2 = EPS * rad * (q_r[2] - q_l[2]);
  dF3 = EPS * rad * (q_r[3] - q_l[3]);

  double half_face = 0.5*mag;

  f[0] = half_face*( (rho_l*V_l             ) + (rho_r*V_r             ) - dF0);
  f[1] = half_face*( (rho_l*u_l*V_l + p_l*r1) + (rho_r*u_r*V_r + p_r*r1) - dF1);
  f[2] = half_face*( (rho_l*v_l*V_l + p_l*r2) + (rho_r*v_r*V_r + p_r*r2) - dF2);
  f[3] = half_face*( (e_l+p_l)*V_l + (e_r+p_r)*V_r                       - dF3);

}


void Adjoint::aflux(){

  int j, k, idx, stride;

//   //
//   // J-direction
//   //
//   for(k=dim->nghost; k< dim->kmax+dim->nghost; k++){
//   for(j=dim->nghost; j<=dim->jmax+dim->nghost; j++){

//     idx = j*dim->jstride + k*dim->kstride;
// #ifdef DO_SIMPLE
//     simpleflux(q[idx-dim->jstride], q[idx], f[idx], dim, grid->Sj[idx]);
// #else
//     roeflux(q[idx-dim->jstride], q[idx], f[idx], dim, grid->Sj[idx]);
// #endif
    

//     rhs[idx-dim->jstride][0] += (f[idx-dim->jstride][0] - f[idx][0])*(j>dim->nghost);
//     rhs[idx-dim->jstride][1] += (f[idx-dim->jstride][1] - f[idx][1])*(j>dim->nghost);
//     rhs[idx-dim->jstride][2] += (f[idx-dim->jstride][2] - f[idx][2])*(j>dim->nghost);
//     rhs[idx-dim->jstride][3] += (f[idx-dim->jstride][3] - f[idx][3])*(j>dim->nghost);

//   }
//   }

//   //
//   // K-direction
//   //
//   for(j=dim->nghost; j< dim->jmax+dim->nghost; j++){
//   for(k=dim->nghost; k<=dim->kmax+dim->nghost; k++){

//     idx = j*dim->jstride + k*dim->kstride;
// #ifdef DO_SIMPLE    
//     simpleflux(q[idx-dim->kstride], q[idx], f[idx], dim, grid->Sk[idx]);
// #else    
//     roeflux(q[idx-dim->kstride], q[idx], f[idx], dim, grid->Sk[idx]);
// #endif    

//     rhs[idx-dim->kstride][0] += (f[idx-dim->kstride][0] - f[idx][0])*(k>dim->nghost);
//     rhs[idx-dim->kstride][1] += (f[idx-dim->kstride][1] - f[idx][1])*(k>dim->nghost);
//     rhs[idx-dim->kstride][2] += (f[idx-dim->kstride][2] - f[idx][2])*(k>dim->nghost);
//     rhs[idx-dim->kstride][3] += (f[idx-dim->kstride][3] - f[idx][3])*(k>dim->nghost);

//   }
//   }

}
