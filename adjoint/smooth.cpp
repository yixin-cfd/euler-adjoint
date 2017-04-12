#include "adjoint.hpp"

//
// Solve a scalar tridiagonal system of equations using the Thomas algorithm
//
//     | D U   |
// A = | L D U |
//     |   L D |
void solve_tridiagonal(double* L, double* D, double* U, double* R,
		       int imax, int istride, int idx_start){
  int idx, idx_m1, idx_p1;
  double inv_main_diag;
  // Initialize for forward sweep
  U[idx_start] = U[idx_start]/D[idx_start];
  R[idx_start] = R[idx_start]/D[idx_start];

  // Forward sweep
  for (int i = 1; i <= imax - 1; i++){
    idx           = idx_start + i*istride;
    idx_m1        = idx - istride;
    inv_main_diag = 1.0/(D[idx] - L[idx]*U[idx_m1]);
    U[idx]        = U[idx]*inv_main_diag;
    R[idx]        = (R[idx] - L[idx]*R[idx_m1])*inv_main_diag;
  }
  // Backward sweep
  for (int i = imax - 2; i >= 0; i--){
    idx    = idx_start + i*istride;
    idx_p1 = idx + istride;
    R[idx] = R[idx] - U[idx]*R[idx_p1];
  }
}

template<int dir>
void fill_LDU(double (*q)[4], double (*L)[4], double (*D)[4], double (*U)[4], 
	      Dim *dim, Grid *grid){

  int i, j, k, idx, stride;
  int jstride = dim->jstride;
  int kstride = dim->kstride;

  double u, v, c;
  double eps, CFL_ratio, spec;
  double (*S)[2];
  double Sx, Sy;

  CFL_ratio = 0.0;

  if(dir == 0){ 
    stride = jstride;
    S = grid->Sj;
  }
  if(dir == 1){ 
    stride = kstride;
    S = grid->Sk;
  }

  memset(L, 0, dim->pts*4*sizeof(double));
  memset(D, 0, dim->pts*4*sizeof(double));
  memset(U, 0, dim->pts*4*sizeof(double));
  // memset(A, 0, dim->pts*4*sizeof(double));

  // int count=0;

  for(k=1; k<dim->ktot-1; k++){
    for(j=1; j<dim->jtot-1; j++){

      idx  = j*jstride + k*kstride;
	   
      Sx   = S[idx][0];
      Sy   = S[idx][1];
	   
      u    = q[idx][1]/q[idx][0];
      v    = q[idx][2]/q[idx][0];
	   
      //     sqrt ( gamma * p / rho ) * ||S||
      c    = sqrt( GAMMA*(GAMMA-1.0)*(q[idx][3]/q[idx][0] - 0.5*(u*u + v*v)) * (Sx*Sx + Sy*Sy) );
      spec = std::abs( Sx*u + Sy*v ) + c;

      eps  = std::max( 0.25*( (CFL_ratio*spec)*(CFL_ratio*spec)-1.0 ), 0.0);

      // if(eps < 1.0e-18){
      // 	count++;
      // }

      for(i=0; i<4; i++){
	L[idx+stride][i] = -eps;
	D[idx       ][i] = 1.0 + 2.0*eps;
	U[idx-stride][i] = -eps;
      }

    }
  }  

  //printf("%d zeros of %d\n", count, (dim->jtot-2)*(dim->ktot-2));

  if(dir == 0){
    // lines in j-dir so loop through k
    for(k=0; k<dim->ktot; k++){
      for(i =0; i<4; i++){
	// first ghost
	L[k*kstride][i]                         = 0.0;
	D[k*kstride][i]                         = 1.0;
	U[k*kstride][i]                         = 0.0;
	// first used L and last used U	          
	L[k*kstride + 1*jstride][i]             = 0.0;
	U[k*kstride + (dim->jtot-2)*jstride][i] = 0.0;
	// last ghost
	L[k*kstride + (dim->jtot-1)*jstride][i] = 0.0;
	D[k*kstride + (dim->jtot-1)*jstride][i] = 1.0;
	U[k*kstride + (dim->jtot-1)*jstride][i] = 0.0;
      }
    }
  }
  if(dir == 1){
    // lines in k-dir so loop through j
    for(j=0; j<dim->jtot; j++){
      for(i =0; i<4; i++){
	// first ghost
	L[j*jstride][i]                         = 0.0;
	D[j*jstride][i]                         = 1.0;
	U[j*jstride][i]                         = 0.0;
	// first used L and last used U	          
	L[j*jstride + 1*kstride][i]             = 0.0;
	U[j*jstride + (dim->ktot-2)*kstride][i] = 0.0;
	// last ghost
	L[j*jstride + (dim->ktot-1)*kstride][i] = 0.0;
	D[j*jstride + (dim->ktot-1)*kstride][i] = 1.0;
	U[j*jstride + (dim->ktot-1)*kstride][i] = 0.0;
      }
    }
  }
  
  
}


void Adjoint::smooth(){

  // double (*saved_rhs)[4] = (double(*)[4])this->scratch;
  // memcpy((void*)saved_rhs, (void*)this->rhs, 4*dim->pts*sizeof(double));

  int start_idx, j, k, var;
  int jstride = dim->jstride;
  int kstride = dim->kstride;

  //
  // J-direction
  //
  fill_LDU<0>(q,L,D,U,dim,grid);

  for(k=1; k<dim->ktot-1; k++){
    for(var=0; var<4; var++){
      start_idx = (1*jstride + k*kstride)*4 + var;
      solve_tridiagonal((double*)L, (double*)D, (double*)U, (double*)rhs, 
			dim->jtot-2, 4*jstride, start_idx);
    }
  }

  //
  // K-direction
  //
  fill_LDU<1>(q,L,D,U,dim,grid);

  for(j=1; j<dim->jtot-1; j++){
    for(var=0; var<4; var++){
      start_idx = (j*jstride + 1*kstride)*4 + var;
      solve_tridiagonal((double*)L, (double*)D, (double*)U, (double*)rhs, 
			dim->ktot-2, 4*kstride, start_idx);
    }
  }


}
