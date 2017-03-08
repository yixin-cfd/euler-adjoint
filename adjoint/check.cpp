#include "adjoint.hpp"
#include "tmpfluxb.h"

void dRdxy(double (*q)[4], double (*psi)[4], double (*xyb)[2], double (*f)[4], 
	   Grid *grid, Dim *dim){

  int j, k, idx, stride;
  double Sxb, Syb;
  double rho, u, v, p, e;
  double dpsi0, dpsi1, dpsi2, dpsi3;
  int mini, maxi;
  double A[4][2];

  //
  // J-direction
  //
  stride = dim->jstride;
  for(k=dim->nghost; k< dim->kmax + dim->nghost; k++){
  for(j=dim->nghost; j<=dim->jmax + dim->nghost; j++){

    idx = j*dim->jstride + k*dim->kstride;

    mini = (j==dim->nghost);
    maxi = (j==dim->nghost + dim->jmax);

    dpsi0 = 0.5*(psi[idx-stride][0]*(1-mini) - psi[idx][0]*(1-maxi));
    dpsi1 = 0.5*(psi[idx-stride][1]*(1-mini) - psi[idx][1]*(1-maxi));
    dpsi2 = 0.5*(psi[idx-stride][2]*(1-mini) - psi[idx][2]*(1-maxi));
    dpsi3 = 0.5*(psi[idx-stride][3]*(1-mini) - psi[idx][3]*(1-maxi));

    // ------------------------------------------------------------------
    // j-1 term
    rho = q[idx-stride][0];
    u   = q[idx-stride][1]/rho;
    v   = q[idx-stride][2]/rho;
    e   = q[idx-stride][3];
    p   = (GAMMA-1.0)*(e - 0.5*rho*(u*u + v*v));

    A[0][0] = rho*u;
    A[1][0] = rho*u*u + p;
    A[2][0] = rho*u*v;
    A[3][0] = (e+p)*u;
    A[0][1] = rho*v;
    A[1][1] = rho*u*v;
    A[2][1] = rho*v*v + p;
    A[3][1] = (e+p)*v;

    // ------------------------------------------------------------------
    // j term
    rho = q[idx][0];
    u   = q[idx][1]/rho;
    v   = q[idx][2]/rho;
    e   = q[idx][3];
    p   = (GAMMA-1.0)*(e - 0.5*rho*(u*u + v*v));

    A[0][0] += rho*u;
    A[1][0] += rho*u*u + p;
    A[2][0] += rho*u*v;
    A[3][0] += (e+p)*u;
    A[0][1] += rho*v;
    A[1][1] += rho*u*v;
    A[2][1] += rho*v*v + p;
    A[3][1] += (e+p)*v;
    // 
    
    Sxb = A[0][0]*dpsi0 + A[1][0]*dpsi1 + A[2][0]*dpsi2 + A[3][0]*dpsi3;
    Syb = A[0][1]*dpsi0 + A[1][1]*dpsi1 + A[2][1]*dpsi2 + A[3][1]*dpsi3;

    // Sj[0] = y[k+1] - y[k  ]
    // Sj[1] = x[k  ] - x[k+1]
    xyb[idx             ][1] -= Sxb;
    xyb[idx+dim->kstride][1] += Sxb;
    xyb[idx             ][0] += Syb;
    xyb[idx+dim->kstride][0] -= Syb;

    // if(j == 5 && k == 1){
    //   // printf("__ %20.14e %20.14e %20.14e %20.14e \n", dpsi0, dpsi1, dpsi2, dpsi3);
    //   // printf("__ %20.14e %20.14e %20.14e %20.14e \n", A[0][0], A[1][0], A[2][0], A[3][0]);
    //   // printf("__ %20.14e \n", A[0][0]*dpsi0 + A[1][0]*dpsi1 + A[2][0]*dpsi2 + A[3][0]*dpsi3);
    //   printf("%d %d: %20.14e %20.14e \n", j, k, xyb[idx][0], xyb[idx][1]);
    // }

  }
  }

  //
  // K-direction
  //
  stride = dim->kstride;
  for(j=dim->nghost; j< dim->jmax + dim->nghost; j++){
  for(k=dim->nghost; k<=dim->kmax + dim->nghost; k++){

    idx = j*dim->jstride + k*dim->kstride;

    mini = (k==dim->nghost);
    maxi = (k==dim->nghost + dim->kmax);

    dpsi0 = 0.5*(psi[idx-stride][0]*(1-mini) - psi[idx][0]*(1-maxi));
    dpsi1 = 0.5*(psi[idx-stride][1]*(1-mini) - psi[idx][1]*(1-maxi));
    dpsi2 = 0.5*(psi[idx-stride][2]*(1-mini) - psi[idx][2]*(1-maxi));
    dpsi3 = 0.5*(psi[idx-stride][3]*(1-mini) - psi[idx][3]*(1-maxi));

    // ------------------------------------------------------------------
    // k-1 term
    rho = q[idx-stride][0];
    u   = q[idx-stride][1]/rho;
    v   = q[idx-stride][2]/rho;
    e   = q[idx-stride][3];
    p   = (GAMMA-1.0)*(e - 0.5*rho*(u*u + v*v));

    A[0][0] = rho*u;
    A[1][0] = rho*u*u + p;
    A[2][0] = rho*u*v;
    A[3][0] = (e+p)*u;
    A[0][1] = rho*v;
    A[1][1] = rho*u*v;
    A[2][1] = rho*v*v + p;
    A[3][1] = (e+p)*v;

    // ------------------------------------------------------------------
    // k term
    rho = q[idx][0];
    u   = q[idx][1]/rho;
    v   = q[idx][2]/rho;
    e   = q[idx][3];
    p   = (GAMMA-1.0)*(e - 0.5*rho*(u*u + v*v));

    A[0][0] += rho*u;
    A[1][0] += rho*u*u + p;
    A[2][0] += rho*u*v;
    A[3][0] += (e+p)*u;
    A[0][1] += rho*v;
    A[1][1] += rho*u*v;
    A[2][1] += rho*v*v + p;
    A[3][1] += (e+p)*v;
    // 
    
    Sxb = A[0][0]*dpsi0 + A[1][0]*dpsi1 + A[2][0]*dpsi2 + A[3][0]*dpsi3;
    Syb = A[0][1]*dpsi0 + A[1][1]*dpsi1 + A[2][1]*dpsi2 + A[3][1]*dpsi3;

    // Sk[0] = y[k  ] - y[k+1]
    // Sk[1] = x[k+1] - x[k  ]
    xyb[idx             ][1] += Sxb;
    xyb[idx+dim->jstride][1] -= Sxb;
    xyb[idx             ][0] -= Syb;
    xyb[idx+dim->jstride][0] += Syb;

  }
  }

}

double Adjoint::check(){

  int i, j, k, idx, idx1;
  double (*f)[4]  = (double (*)[4])this->scratch;

  double f1[2], f2[2], d1[4], d2[4];
  int mini, maxi;

  // printf("\nmanually changing things ... \n\n");
  // for(k=0; k<dim->ktot; k++){
  //   for(j=0; j<dim->jtot; j++){
  //     idx  = j*dim->jstride + k*dim->kstride;
  //     // if((j == 6 || j == 5) && k == 1){
  //     if(j == 5 && k == 1){
  // 	psi[idx][0] = 1.0;
  // 	// psi[idx][0] = psi[idx][0];
  //     } else {
  // 	psi[idx][0] = 0.0;
  //     }
  //     psi[idx][1] = 0.0;
  //     psi[idx][2] = 0.0;
  //     psi[idx][3] = 0.0;
  //   }
  // }


  memset(xyb, 0, 2*dim->pts*sizeof(double));
  dRdxy(q, psi, xyb, f, grid, dim);

  // memset(xyb, 0, 2*dim->pts*sizeof(double));
  // for(k=dim->nghost; k< dim->kmax+dim->nghost; k++){
  // for(j=dim->nghost; j<=dim->jmax+dim->nghost; j++){
  //   idx  = j*dim->jstride + k*dim->kstride;
  //   idx1 = idx-dim->jstride;
  //   mini = (j == dim->nghost);
  //   maxi = (j == dim->nghost + dim->jmax);
  //   f1[0] = 0.0; f1[1] = 0.0; f1[2] = 0.0; f1[3] = 0.0;
  //   f2[0] = 0.0; f2[1] = 0.0; f2[2] = 0.0; f2[3] = 0.0;
  //   tmpflux_b(q[idx1], q[idx], d1, psi[idx1], d2, psi[idx], dim, 
  // 	      grid->xy[idx], f1,
  // 	      grid->xy[idx+dim->kstride], f2,
  // 	      mini, maxi);
  //   xyb[idx             ][0] += f1[0];
  //   xyb[idx             ][1] += f1[1];
  //   xyb[idx+dim->kstride][0] += f2[0];
  //   xyb[idx+dim->kstride][1] += f2[1];
  // }
  // }
  // for(j=dim->nghost; j< dim->jmax+dim->nghost; j++){
  // for(k=dim->nghost; k<=dim->kmax+dim->nghost; k++){
  //   idx  = j*dim->jstride + k*dim->kstride;
  //   idx1 = idx-dim->kstride;
  //   mini = (k == dim->nghost);
  //   maxi = (k == dim->nghost + dim->kmax);
  //   f1[0] = 0.0; f1[1] = 0.0; f1[2] = 0.0; f1[3] = 0.0;
  //   f2[0] = 0.0; f2[1] = 0.0; f2[2] = 0.0; f2[3] = 0.0;
  //   tmpflux_b(q[idx1], q[idx], d1, psi[idx1], d2, psi[idx], dim, 
  // 	      grid->xy[idx+dim->jstride], xyb[idx+dim->jstride],
  // 	      grid->xy[idx], xyb[idx],
  // 	      mini, maxi);
  //   xyb[idx1][0] += f1[0];
  //   xyb[idx1][1] += f1[1];
  //   xyb[idx][0]  += f2[0];
  //   xyb[idx][1]  += f2[1];
  // }
  // }


  // double tmpres = 0.0;
  // for(k=0; k<dim->ktot; k++){
  //   for(j=0; j<dim->jtot; j++){
  //     idx = j*dim->jstride + k*dim->kstride;
  //     tmpres += xyb[idx][0]*xyb[idx][0] + xyb[idx][1]*xyb[idx][1];
  //     // if(std::abs(xyb[idx][1]) > 1e-12){
  //     if(j == 5 && k == 1){
  //     	printf("%d %d: %20.14e %20.14e \n", j, k, xyb[idx][0], xyb[idx][1]);
  //     }
  //   }
  // }
  // printf("x res: %20.14e\n", tmpres);

  double (*xd)[2] = (double (*)[2])this->scratch;

  for(i=0; i<dim->pts; i++){
    xd[i][0] =  grid->xy[i][1];
    xd[i][1] = -grid->xy[i][0];
  }

  double liftd = 0.0;
  for(i=0; i<dim->pts; i++){

    liftd += xyb[i][0]*xd[i][0] + xyb[i][1]*xd[i][1];

  }

  return liftd;

}
