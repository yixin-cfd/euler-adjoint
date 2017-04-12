#include "adjoint.hpp"
#include "tmpfluxb.h"

#define EPS 0.25

void dRdxy(double (*q)[4], double (*psi)[4], double (*xyb)[2], double (*f)[4], 
	   Grid *grid, Dim *dim){

  int j, k, idx, stride;
  double Sxb, Syb;
  double rho, u, v, p, e, c, V;
  double dpsi0, dpsi1, dpsi2, dpsi3;
  int mini, maxi;
  double A[4][2];
  double dspec_x, dspec_y, imag;

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

    imag  = 1.0/sqrt(grid->Sj[idx][0]*grid->Sj[idx][0] + grid->Sj[idx][1]*grid->Sj[idx][1]);

    // ------------------------------------------------------------------
    // j-1 term
    rho      = q[idx-stride][0];
    u        = q[idx-stride][1]/rho;
    v        = q[idx-stride][2]/rho;
    e        = q[idx-stride][3];
    p        = (GAMMA-1.0)*(e - 0.5*rho*(u*u + v*v));
    c        = sqrt(GAMMA*p/rho);
    V        = grid->Sj[idx][0]*u + grid->Sj[idx][1]*v;

    if(V >= 0.0){
      dspec_x  = u + grid->Sj[idx][0]*imag*c;
      dspec_y  = v + grid->Sj[idx][1]*imag*c;
    } else {
      dspec_x  = -u + grid->Sj[idx][0]*imag*c;
      dspec_y  = -v + grid->Sj[idx][1]*imag*c;
    }
	     
    A[0][0]  = rho*u;
    A[1][0]  = rho*u*u + p;
    A[2][0]  = rho*u*v;
    A[3][0]  = (e+p)*u;
    A[0][1]  = rho*v;
    A[1][1]  = rho*u*v;
    A[2][1]  = rho*v*v + p;
    A[3][1]  = (e+p)*v;

    // ------------------------------------------------------------------
    // j term
    rho      = q[idx][0];
    u        = q[idx][1]/rho;
    v        = q[idx][2]/rho;
    e        = q[idx][3];
    p        = (GAMMA-1.0)*(e - 0.5*rho*(u*u + v*v));
    c        = sqrt(GAMMA*p/rho);
    V        = grid->Sj[idx][0]*u + grid->Sj[idx][1]*v;

    if(V >= 0.0){
      dspec_x  += u + grid->Sj[idx][0]*imag*c;
      dspec_y  += v + grid->Sj[idx][1]*imag*c;
    } else {
      dspec_x  += -u + grid->Sj[idx][0]*imag*c;
      dspec_y  += -v + grid->Sj[idx][1]*imag*c;
    }

    A[0][0] += rho*u;
    A[1][0] += rho*u*u + p;
    A[2][0] += rho*u*v;
    A[3][0] += (e+p)*u;
    A[0][1] += rho*v;
    A[1][1] += rho*u*v;
    A[2][1] += rho*v*v + p;
    A[3][1] += (e+p)*v;

    // ------------------------------------------------------------------
    // dissipation

    A[0][0] += -0.5*EPS*dspec_x*(q[idx][0]-q[idx-stride][0]);
    A[1][0] += -0.5*EPS*dspec_x*(q[idx][1]-q[idx-stride][1]);
    A[2][0] += -0.5*EPS*dspec_x*(q[idx][2]-q[idx-stride][2]);
    A[3][0] += -0.5*EPS*dspec_x*(q[idx][3]-q[idx-stride][3]);
    A[0][1] += -0.5*EPS*dspec_y*(q[idx][0]-q[idx-stride][0]);
    A[1][1] += -0.5*EPS*dspec_y*(q[idx][1]-q[idx-stride][1]);
    A[2][1] += -0.5*EPS*dspec_y*(q[idx][2]-q[idx-stride][2]);
    A[3][1] += -0.5*EPS*dspec_y*(q[idx][3]-q[idx-stride][3]);

    // ---------------

    Sxb = A[0][0]*dpsi0 + A[1][0]*dpsi1 + A[2][0]*dpsi2 + A[3][0]*dpsi3;
    Syb = A[0][1]*dpsi0 + A[1][1]*dpsi1 + A[2][1]*dpsi2 + A[3][1]*dpsi3;

    // Sj[0] = y[k+1] - y[k  ]
    // Sj[1] = x[k  ] - x[k+1]
    xyb[idx             ][1] -= Sxb;
    xyb[idx+dim->kstride][1] += Sxb;
    xyb[idx             ][0] += Syb;
    xyb[idx+dim->kstride][0] -= Syb;

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

    imag  = 1.0/sqrt(grid->Sk[idx][0]*grid->Sk[idx][0] + grid->Sk[idx][1]*grid->Sk[idx][1]);

    dpsi0 = 0.5*(psi[idx-stride][0]*(1-mini) - psi[idx][0]*(1-maxi));
    dpsi1 = 0.5*(psi[idx-stride][1]*(1-mini) - psi[idx][1]*(1-maxi));
    dpsi2 = 0.5*(psi[idx-stride][2]*(1-mini) - psi[idx][2]*(1-maxi));
    dpsi3 = 0.5*(psi[idx-stride][3]*(1-mini) - psi[idx][3]*(1-maxi));

    // ------------------------------------------------------------------
    // k-1 term
    rho      = q[idx-stride][0];
    u        = q[idx-stride][1]/rho;
    v        = q[idx-stride][2]/rho;
    e        = q[idx-stride][3];
    p        = (GAMMA-1.0)*(e - 0.5*rho*(u*u + v*v));
    c        = sqrt(GAMMA*p/rho);
    V        = grid->Sk[idx][0]*u + grid->Sk[idx][1]*v;

    if(V >= 0.0){
      dspec_x  = u + grid->Sk[idx][0]*imag*c;
      dspec_y  = v + grid->Sk[idx][1]*imag*c;
    } else {
      dspec_x  = -u + grid->Sk[idx][0]*imag*c;
      dspec_y  = -v + grid->Sk[idx][1]*imag*c;
    }

    A[0][0]  = rho*u;
    A[1][0]  = rho*u*u + p;
    A[2][0]  = rho*u*v;
    A[3][0]  = (e+p)*u;
    A[0][1]  = rho*v;
    A[1][1]  = rho*u*v;
    A[2][1]  = rho*v*v + p;
    A[3][1]  = (e+p)*v;

    // ------------------------------------------------------------------
    // k term
    rho      = q[idx][0];
    u        = q[idx][1]/rho;
    v        = q[idx][2]/rho;
    e        = q[idx][3];
    p        = (GAMMA-1.0)*(e - 0.5*rho*(u*u + v*v));
    c        = sqrt(GAMMA*p/rho);
    V        = grid->Sk[idx][0]*u + grid->Sk[idx][1]*v;

    if(V >= 0.0){
      dspec_x  += u + grid->Sk[idx][0]*imag*c;
      dspec_y  += v + grid->Sk[idx][1]*imag*c;
    } else {   
      dspec_x  += -u + grid->Sk[idx][0]*imag*c;
      dspec_y  += -v + grid->Sk[idx][1]*imag*c;
    }

    A[0][0] += rho*u;
    A[1][0] += rho*u*u + p;
    A[2][0] += rho*u*v;
    A[3][0] += (e+p)*u;
    A[0][1] += rho*v;
    A[1][1] += rho*u*v;
    A[2][1] += rho*v*v + p;
    A[3][1] += (e+p)*v;

    // ------------------------------------------------------------------
    // dissipation

    A[0][0] += -0.5*EPS*dspec_x*(q[idx][0]-q[idx-stride][0]);
    A[1][0] += -0.5*EPS*dspec_x*(q[idx][1]-q[idx-stride][1]);
    A[2][0] += -0.5*EPS*dspec_x*(q[idx][2]-q[idx-stride][2]);
    A[3][0] += -0.5*EPS*dspec_x*(q[idx][3]-q[idx-stride][3]);
    A[0][1] += -0.5*EPS*dspec_y*(q[idx][0]-q[idx-stride][0]);
    A[1][1] += -0.5*EPS*dspec_y*(q[idx][1]-q[idx-stride][1]);
    A[2][1] += -0.5*EPS*dspec_y*(q[idx][2]-q[idx-stride][2]);
    A[3][1] += -0.5*EPS*dspec_y*(q[idx][3]-q[idx-stride][3]);


    // --------
    
    Sxb = A[0][0]*dpsi0 + A[1][0]*dpsi1 + A[2][0]*dpsi2 + A[3][0]*dpsi3;
    Syb = A[0][1]*dpsi0 + A[1][1]*dpsi1 + A[2][1]*dpsi2 + A[3][1]*dpsi3;


    xyb[idx             ][1] += Sxb;
    xyb[idx+dim->jstride][1] -= Sxb;
    xyb[idx             ][0] -= Syb;
    xyb[idx+dim->jstride][0] += Syb;

  }
  }
}

void dBCdxy(BC bc, double (*q)[4], double (*rhs)[4], double (*xyb)[2], Dim *dim, Grid *grid){

  int j, k, idx, midx, widx;
  double Sx, Sy, Sxb, Syb, SS;
  int jstride, kstride;
  jstride = dim->jstride;
  kstride = dim->kstride;

  for(k=bc.ks; k<=bc.ke; k++){
    for(j=bc.js; j<=bc.je; j++){
	
      idx  = j*dim->jstride + k*dim->kstride;
      midx = j*jstride + (2*dim->nghost - k - 1)*kstride;
      widx = j*jstride + dim->nghost*kstride;
  
      Sx   = grid->Sk[widx][0];
      Sy   = grid->Sk[widx][1];
	   
      SS   = Sx*Sx + Sy*Sy;

      Sxb  = ( -4.0*q[midx][1]*( (SS * Sx -     Sx*Sx*Sx)/(SS*SS))
      	       -2.0*q[midx][2]*( (SS * Sy - 2.0*Sx*Sx*Sy)/(SS*SS)) )*rhs[idx][1];
      Syb  = ( -4.0*q[midx][1]*( (        -     Sx*Sx*Sy)/(SS*SS))
      	       -2.0*q[midx][2]*( (SS * Sx - 2.0*Sx*Sy*Sy)/(SS*SS)) )*rhs[idx][1];

      Sxb += ( -2.0*q[midx][1]*( (SS * Sy - 2.0*Sx*Sx*Sy)/(SS*SS))
      	       -4.0*q[midx][2]*( (        -     Sy*Sy*Sx)/(SS*SS)) )*rhs[idx][2];
      Syb += ( -2.0*q[midx][1]*( (SS * Sx - 2.0*Sx*Sy*Sy)/(SS*SS))
      	       -4.0*q[midx][2]*( (SS * Sy -     Sy*Sy*Sy)/(SS*SS)) )*rhs[idx][2];

      xyb[midx             ][1] += Sxb;
      xyb[midx+dim->jstride][1] -= Sxb;
      xyb[midx             ][0] -= Syb;
      xyb[midx+dim->jstride][0] += Syb;

    }
  }


}

double Adjoint::check(){

  int i, j, k, idx, idx1;
  double (*f)[4]  = (double (*)[4])this->scratch;
  double f1[2], f2[2], d1[2], d2[2];
  int jstride = dim->jstride;
  int kstride = dim->kstride;

  int mini, maxi;

  memset( rhs, 0, 4*dim->pts*sizeof(double));
  memset( xyb, 0, 2*dim->pts*sizeof(double));

  // fill rhs values WITHOUT dependence on cost function
  this->aflux();

  // find dependence of residual on x
  dRdxy(q, psi, xyb, f, grid, dim);
  
  // find dependence of bcs on x
  for(i=0; i<euler->inputs->nbc; i++){
    if (euler->bc[i].type != WALL_BC) continue;
    dBCdxy(euler->bc[i], q, rhs, xyb, dim, grid);
  }

  // normally here we would multiply xyb by xyd ( a delta grid change
  // ). Since we dont have that, we can make one up so that we can
  // still get a feel for the convergence of the sensitivity.
  double liftd = 0.0;
  double tmp1, tmp2;
  for(i=0; i<dim->pts; i++){
    tmp1 = grid->xy[i][1];
    tmp2 = -grid->xy[i][0];
    liftd += xyb[i][0]*tmp1 + xyb[i][1]*tmp2;
  }

  return liftd;

}
