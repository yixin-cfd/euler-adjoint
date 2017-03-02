#include "adjoint.hpp"
#include "tmpfluxb.h"

// #define DO_SIMPLE
#define EPS 0.25

void adj_flux(double* q_l, double* q_r, double* psi_l, double* psi_r, double* f, 
	      Dim *dim, double S[2], double A[4][4]){

  double Sx, Sy, phi, u, v, V, E;
  double a1;
  double a2 = GAMMA-1;
  double a3 = GAMMA-2;
  double dpsi0, dpsi1, dpsi2, dpsi3;

  Sx         = S[0];
  Sy         = S[1];


  dpsi0      = psi_r[0] - psi_l[0];
  dpsi1      = psi_r[1] - psi_l[1];
  dpsi2      = psi_r[2] - psi_l[2];
  dpsi3      = psi_r[3] - psi_l[3];

  u          = 0.5*(q_l[1]/q_l[0] + q_r[1]/q_r[0]);
  v          = 0.5*(q_l[2]/q_l[0] + q_r[2]/q_r[0]);
  E          = 0.5*(q_l[3]/q_l[0] + q_r[3]/q_r[0]);
  
  // // ---------------------------------------------------------------------
  // // LEFT
  // //
  // u          = q_l[1]/q_l[0];
  // v          = q_l[2]/q_l[0];
  // E          = q_l[3]/q_l[0];

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
  f[0]  = A[0][0]*dpsi0 + A[1][0]*dpsi1 + A[2][0]*dpsi2 + A[3][0]*dpsi3;
  f[1]  = A[0][1]*dpsi0 + A[1][1]*dpsi1 + A[2][1]*dpsi2 + A[3][1]*dpsi3;
  f[2]  = A[0][2]*dpsi0 + A[1][2]*dpsi1 + A[2][2]*dpsi2 + A[3][2]*dpsi3;
  f[3]  = A[0][3]*dpsi0 + A[1][3]*dpsi1 + A[2][3]*dpsi2 + A[3][3]*dpsi3;
  // // ---------------------------------------------------------------------
  // // RIGHT
  // //
  // u          = q_r[1]/q_r[0];
  // v          = q_r[2]/q_r[0];
  // E          = q_r[3]/q_r[0];

  // V          = Sx*u + Sy*v;
  // c          = sqrt( a2*GAMMA*(E - 0.5*(u*u + v*v)) * (Sx*Sx + Sy*Sy) );
  // phi        = 0.5*a2*(u*u + v*v);
  // a1         = GAMMA*E - phi;
  // // build the Jacobian
  // A[0][0]    = 0.0;
  // A[0][1]    = Sx;
  // A[0][2]    = Sy;
  // A[0][3]    = 0.0;
  // A[1][0]    = Sx*phi - u*V;
  // A[1][1]    = V-a3*Sx*u;
  // A[1][2]    = Sy*u-a2*Sx*v;
  // A[1][3]    = a2*Sx;
  // A[2][0]    = Sy*phi-v*V;
  // A[2][1]    = Sx*v-a2*Sy*u;
  // A[2][2]    = V-a3*Sy*v;
  // A[2][3]    = a2*Sy;
  // A[3][0]    = V*(phi-a1);
  // A[3][1]    = Sx*a1-a2*u*V;
  // A[3][2]    = Sy*a1-a2*v*V;
  // A[3][3]    = GAMMA*V;
  // // // // We want  [ A^T ]*[dpsi] - dissipation
  // f[0] += A[0][0]*dpsi0 + A[1][0]*dpsi1 + A[2][0]*dpsi2 + A[3][0]*dpsi3;
  // f[1] += A[0][1]*dpsi0 + A[1][1]*dpsi1 + A[2][1]*dpsi2 + A[3][1]*dpsi3;
  // f[2] += A[0][2]*dpsi0 + A[1][2]*dpsi1 + A[2][2]*dpsi2 + A[3][2]*dpsi3;
  // f[3] += A[0][3]*dpsi0 + A[1][3]*dpsi1 + A[2][3]*dpsi2 + A[3][3]*dpsi3;

  // f[0]  = A[0][0]*dpsi0 + A[0][1]*dpsi0 + A[0][2]*dpsi0 + A[0][3]*dpsi0;
  // f[1]  = A[1][0]*dpsi1 + A[1][1]*dpsi1 + A[1][2]*dpsi1 + A[1][3]*dpsi1;
  // f[2]  = A[2][0]*dpsi2 + A[2][1]*dpsi2 + A[2][2]*dpsi2 + A[2][3]*dpsi2;
  // f[3]  = A[3][0]*dpsi3 + A[3][1]*dpsi3 + A[3][2]*dpsi3 + A[3][3]*dpsi3;

  // f[0]  = A[0][0]*dpsi0 + A[0][1]*dpsi1 + A[0][2]*dpsi2 + A[0][3]*dpsi3;
  // f[1]  = A[1][0]*dpsi0 + A[1][1]*dpsi1 + A[1][2]*dpsi2 + A[1][3]*dpsi3;
  // f[2]  = A[2][0]*dpsi0 + A[2][1]*dpsi1 + A[2][2]*dpsi2 + A[2][3]*dpsi3;
  // f[3]  = A[3][0]*dpsi0 + A[3][1]*dpsi1 + A[3][2]*dpsi2 + A[3][3]*dpsi3;

  // f[0]  = A[0][0]*dpsi0 + A[1][0]*dpsi0 + A[2][0]*dpsi0 + A[3][0]*dpsi0;
  // f[1]  = A[0][1]*dpsi1 + A[1][1]*dpsi1 + A[2][1]*dpsi1 + A[3][1]*dpsi1;
  // f[2]  = A[0][2]*dpsi2 + A[1][2]*dpsi2 + A[2][2]*dpsi2 + A[3][2]*dpsi2;
  // f[3]  = A[0][3]*dpsi3 + A[1][3]*dpsi3 + A[2][3]*dpsi3 + A[3][3]*dpsi3;

  // f[0] *= 0.5;
  // f[1] *= 0.5;
  // f[2] *= 0.5;
  // f[3] *= 0.5;

}

void adj_diss(double* q_l, double* q_r, double* psi_l, double* psi_r, double* d,
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

  int j, k, idx, idx1;
  
  double (*f)[4] = (double (*)[4])this->scratch;

  double (*g)[4] = new double[dim->pts][4];

  double (*d)[4] = (double (*)[4])this->scratch;

  double f1[4], f2[4], d1[4], d2[4];

  double A[4][4];

  double fac1 = 0.5;
  double fac2 = 1.0;

  int mini, maxi;

  // if(step_number == 1){
  //   printf("\nmanually changing things ... \n\n");
  //   for(k=0; k<dim->ktot; k++){
  //     for(j=0; j<dim->jtot; j++){
  // 	idx  = j*dim->jstride + k*dim->kstride;
  // 	psi[idx][0] = 1.0*(j == 5 && k == 1);
  // 	psi[idx][1] = 0.0;
  // 	psi[idx][2] = 0.0;
  // 	psi[idx][3] = 0.0;
  //     }
  //   }
  // }

  memset(g, 0, 4*dim->pts*sizeof(double));

  //
  // J-Direction
  //
  for(k=dim->nghost; k< dim->kmax+dim->nghost; k++){
  for(j=dim->nghost; j<=dim->jmax+dim->nghost; j++){

    idx  = j*dim->jstride + k*dim->kstride;
    idx1 = idx-dim->jstride;

    // adj_flux(q[idx1], q[idx], psi[idx1], psi[idx], f[idx], dim, grid->Sj[idx], A);
    // rhs[idx1][0] += fac1*(f[idx1][0] + f[idx][0])*(j>dim->nghost);
    // rhs[idx1][1] += fac1*(f[idx1][1] + f[idx][1])*(j>dim->nghost);
    // rhs[idx1][2] += fac1*(f[idx1][2] + f[idx][2])*(j>dim->nghost);
    // rhs[idx1][3] += fac1*(f[idx1][3] + f[idx][3])*(j>dim->nghost);
    // rhs[idx][0] += fac1*(f[idx][0])*(j>dim->nghost);
    // rhs[idx][1] += fac1*(f[idx][1])*(j>dim->nghost);
    // rhs[idx][2] += fac1*(f[idx][2])*(j>dim->nghost);
    // rhs[idx][3] += fac1*(f[idx][3])*(j>dim->nghost);

    // --------------------------------------------------------------------------------
    
    f1[0] = 0.0; f1[1] = 0.0; f1[2] = 0.0; f1[3] = 0.0;
    f2[0] = 0.0; f2[1] = 0.0; f2[2] = 0.0; f2[3] = 0.0;
    mini = (j == dim->nghost);
    maxi = (j == dim->nghost + dim->jmax);
    tmpflux_b(q[idx1], f1, q[idx], f2, d1, psi[idx1], d2, psi[idx], 
    	      dim, grid->Sj[idx], mini, maxi);
    rhs[idx1][0] += f1[0];
    rhs[idx1][1] += f1[1];
    rhs[idx1][2] += f1[2];
    rhs[idx1][3] += f1[3];
    rhs[idx][0]  += f2[0];
    rhs[idx][1]  += f2[1];
    rhs[idx][2]  += f2[2];
    rhs[idx][3]  += f2[3];

    // --------------------------------------------------------------------------------

    // tmpflux_b(q[idx1], rhs[idx1], q[idx], rhs[idx], d1, psi[idx1], d2, psi[idx], 
    // 	      dim, grid->Sj[idx], mini, maxi);


    // tmpflux_b(q[idx1], f1, q[idx], f2, d1, d2, dim, grid->Sj[idx]);


    // if(j == 5 && k == 1){
    //   printf("__ %20.16e %20.16e %20.16e %20.16e \n", f[idx1][0],f[idx1][1],f[idx1][2],f[idx1][3]);
    //   printf("__ %20.16e %20.16e %20.16e %20.16e \n",  f[idx][0], f[idx][1], f[idx][2], f[idx][3]);
    //   printf("__ %20.16e %20.16e %20.16e %20.16e \n",     f1[0],     f1[1],     f1[2],     f1[3]);
    //   printf("__ %20.16e %20.16e %20.16e %20.16e \n",     f2[0],     f2[1],     f2[2],     f2[3]);
    //   printf("------- ( %20.16e, %20.16e ) ------\n", grid->Sj[idx][0], grid->Sj[idx][1]);
    // }

  }
  }

  //
  // K-direction
  //
  for(j=dim->nghost; j< dim->jmax+dim->nghost; j++){
  for(k=dim->nghost; k<=dim->kmax+dim->nghost; k++){

    idx  = j*dim->jstride + k*dim->kstride;
    idx1 = idx-dim->kstride;

    // adj_flux(q[idx1], q[idx], psi[idx1], psi[idx], f[idx], dim, grid->Sk[idx], A);
    // rhs[idx1][0] += fac1*(f[idx1][0] + f[idx][0])*(k>dim->nghost);
    // rhs[idx1][1] += fac1*(f[idx1][1] + f[idx][1])*(k>dim->nghost);
    // rhs[idx1][2] += fac1*(f[idx1][2] + f[idx][2])*(k>dim->nghost);
    // rhs[idx1][3] += fac1*(f[idx1][3] + f[idx][3])*(k>dim->nghost);
    // rhs[idx][0] += fac1*(f[idx][0])*(j>dim->nghost);
    // rhs[idx][1] += fac1*(f[idx][1])*(j>dim->nghost);
    // rhs[idx][2] += fac1*(f[idx][2])*(j>dim->nghost);
    // rhs[idx][3] += fac1*(f[idx][3])*(j>dim->nghost);

    // ----------------------------------------------------------------------------------------

    f1[0] = 0.0; f1[1] = 0.0; f1[2] = 0.0; f1[3] = 0.0;
    f2[0] = 0.0; f2[1] = 0.0; f2[2] = 0.0; f2[3] = 0.0;
    mini = (k == dim->nghost);
    maxi = (k == dim->nghost + dim->kmax);
    tmpflux_b(q[idx1], f1, q[idx], f2, d1, psi[idx1], d2, psi[idx], 
    	      dim, grid->Sk[idx], mini, maxi);
    rhs[idx1][0] += f1[0];
    rhs[idx1][1] += f1[1];
    rhs[idx1][2] += f1[2];
    rhs[idx1][3] += f1[3];
    rhs[idx][0]  += f2[0];
    rhs[idx][1]  += f2[1];
    rhs[idx][2]  += f2[2];
    rhs[idx][3]  += f2[3];

    // ----------------------------------------------------------------------------------------

    // tmpflux_b(q[idx1], rhs[idx1], q[idx], rhs[idx], d1, psi[idx1], d2, psi[idx], 
    // 	      dim, grid->Sk[idx], mini, maxi);
    // tmpflux_b(q[idx1], f1, q[idx], f2, d1, d2, dim, grid->Sk[idx]);


  }
  }

  // //
  // // Dissipation: J-Direction 
  // //
  // for(k=dim->nghost; k< dim->kmax+dim->nghost; k++){
  // for(j=dim->nghost; j<=dim->jmax+dim->nghost; j++){

  //   idx  = j*dim->jstride + k*dim->kstride;
  //   idx1 = idx-dim->jstride;

  //   adj_diss(q[idx1], q[idx], psi[idx1], psi[idx], d[idx], dim, grid->Sj[idx]);
  //   rhs[idx1][0] += fac2*(d[idx1][0] - d[idx][0])*(j>dim->nghost);
  //   rhs[idx1][1] += fac2*(d[idx1][1] - d[idx][1])*(j>dim->nghost);
  //   rhs[idx1][2] += fac2*(d[idx1][2] - d[idx][2])*(j>dim->nghost);
  //   rhs[idx1][3] += fac2*(d[idx1][3] - d[idx][3])*(j>dim->nghost);
  // }
  // }
  // //
  // // Dissipation: K-direction
  // //
  // for(j=dim->nghost; j< dim->jmax+dim->nghost; j++){
  // for(k=dim->nghost; k<=dim->kmax+dim->nghost; k++){

  //   idx  = j*dim->jstride + k*dim->kstride;
  //   idx1 = idx-dim->kstride;

  //   adj_diss(q[idx1], q[idx], psi[idx1], psi[idx], d[idx], dim, grid->Sk[idx]);
  //   rhs[idx1][0] += fac2*(d[idx1][0] - d[idx][0])*(k>dim->nghost);
  //   rhs[idx1][1] += fac2*(d[idx1][1] - d[idx][1])*(k>dim->nghost);
  //   rhs[idx1][2] += fac2*(d[idx1][2] - d[idx][2])*(k>dim->nghost);
  //   rhs[idx1][3] += fac2*(d[idx1][3] - d[idx][3])*(k>dim->nghost);
  // }
  // }


  delete g;

}
