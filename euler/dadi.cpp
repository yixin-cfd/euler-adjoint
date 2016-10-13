#include "euler.hpp"

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

void compute_T_inv(double q[4], double nx, double ny, double T_inv[16]){

  double rho, u, v, p, c, irho, phi_sq, alpha, beta, theta, ic2;
  rho    = q[0];
  irho   = 1.0 / rho;
  u      = q[1]*irho;
  v      = q[2]*irho;

  phi_sq = 0.5*(GAMMA-1.0)*(u*u + v*v);

  p      = (GAMMA-1.0)*q[3] - rho*phi_sq;
  c      = sqrt(GAMMA*p*irho);
  ic2    = 1.0 / (c*c);

  alpha  = rho / (sqrt(2.0)*c);
  beta   = alpha * irho * irho;
  theta  = nx*u + ny*v;

  T_inv[ 0] = 1.0 - phi_sq*ic2;
  T_inv[ 1] = (GAMMA-1.0)*ic2*u;
  T_inv[ 2] = (GAMMA-1.0)*ic2*v;
  T_inv[ 3] = -(GAMMA-1.0)*ic2;

  T_inv[ 4] = -irho*(ny*u - nx*v);
  T_inv[ 5] = ny*irho;
  T_inv[ 6] = -nx*irho;
  T_inv[ 7] = 0;

  T_inv[ 8] = beta*(phi_sq - c*theta);
  T_inv[ 9] = beta*(nx*c-(GAMMA-1.0)*u);
  T_inv[10] = beta*(ny*c-(GAMMA-1.0)*v);
  T_inv[11] = beta*(GAMMA-1.0);

  T_inv[12] = beta*(phi_sq + c*theta);
  T_inv[13] = -beta*(nx*c+(GAMMA-1.0)*u);
  T_inv[14] = -beta*(ny*c+(GAMMA-1.0)*v);
  T_inv[15] = beta*(GAMMA-1.0);
  
}

void compute_T(double q[4], double nx, double ny, double T[16]){
  double rho, u, v, p, c, irho, phi_sq, alpha, beta, theta, ic2;
  rho    = q[0];
  irho   = 1.0 / rho;
  u      = q[1]*irho;
  v      = q[2]*irho;

  phi_sq = 0.5*(GAMMA-1.0)*(u*u + v*v);

  p      = (GAMMA-1.0)*q[3] - rho*phi_sq;
  c      = sqrt(GAMMA*p*irho);
  ic2    = 1.0 / (c*c);

  alpha  = rho / (sqrt(2.0)*c);
  // beta   = alpha * irho * irho;
  theta  = nx*u + ny*v;

  T[ 0] = 1.0;
  T[ 1] = 0;
  T[ 2] = alpha;
  T[ 3] = alpha;

  T[ 4] = u;		 
  T[ 5] = ny*rho;	 
  T[ 6] = alpha*(u+nx*c);
  T[ 7] = alpha*(u-nx*c);

  T[ 8] = v;		 
  T[ 9] = -nx*rho;	 
  T[10] = alpha*(v+ny*c);
  T[11] = alpha*(v-ny*c);

  T[12] = phi_sq / (GAMMA-1.0);
  T[13] = rho*(ny*u-nx*v);
  T[14] = alpha*( (phi_sq + c*c) / (GAMMA-1.0) + c*theta);
  T[15] = alpha*( (phi_sq + c*c) / (GAMMA-1.0) - c*theta);

}

void Tk_inv_Tl(double kx, double ky, double lx, double ly, double N[16]){

  double u = 1.0 / sqrt(2.0);
  double m1 = kx*lx + ky*ly;
  double m2 = kx*ly - ky*lx;

  N[ 0] = 1.0;
  N[ 1] = 0.0;
  N[ 2] = 0.0;
  N[ 3] = 0.0;

  N[ 4] = 0.0;
  N[ 5] = m1;
  N[ 6] = -u*m2;
  N[ 7] = u*m2;

  N[ 8] = 0.0;
  N[ 9] = u*m2;
  N[10] = u*u*(1+m1);
  N[11] = u*u*(1-m1);

  N[12] = 0.0;
  N[13] = -u*m2;
  N[14] = u*u*(1-m1);
  N[15] = u*u*(1+m1);
  
}


void invert_xi(double (*rhs)[4], double (*q)[4], Dim *dim, Grid *grid){

  int i, j, k;
  double nx, ny, tmp;
  double T_inv[16];
  double r0, r1, r2, r3;

  for(k=dim->nghost; k<dim->kmax + dim->nghost; k++){
  for(j=dim->nghost; j<dim->jmax + dim->nghost; j++){

    i = j*dim->jstride + k*dim->kstride;

    nx = grid->Sj[i][0];
    ny = grid->Sj[i][1];
    tmp = 1.0 / sqrt(nx*nx + ny*ny);
    nx = nx * tmp;
    ny = ny * tmp;

    compute_T_inv(q[i], nx, ny, T_inv);

    r0 = T_inv[ 0]*rhs[i][0] + T_inv[ 1]*rhs[i][1] + T_inv[ 2]*rhs[i][2] + T_inv[ 3]*rhs[i][3];
    r1 = T_inv[ 4]*rhs[i][0] + T_inv[ 5]*rhs[i][1] + T_inv[ 6]*rhs[i][2] + T_inv[ 7]*rhs[i][3];
    r2 = T_inv[ 8]*rhs[i][0] + T_inv[ 9]*rhs[i][1] + T_inv[10]*rhs[i][2] + T_inv[11]*rhs[i][3];
    r3 = T_inv[12]*rhs[i][0] + T_inv[13]*rhs[i][1] + T_inv[14]*rhs[i][2] + T_inv[15]*rhs[i][3];

    rhs[i][0] = r0;
    rhs[i][1] = r1;
    rhs[i][2] = r2;
    rhs[i][3] = r3;
  }
  }
  
}

void invert_eta(double (*rhs)[4], double (*q)[4], Dim *dim, Grid *grid){

  int i, j, k;
  double nx, ny, tmp;
  double T[16];
  double r0, r1, r2, r3;

  for(k=dim->nghost; k<dim->kmax + dim->nghost; k++){
  for(j=dim->nghost; j<dim->jmax + dim->nghost; j++){

    i = j*dim->jstride + k*dim->kstride;

    nx = grid->Sk[i][0];
    ny = grid->Sk[i][1];
    tmp = 1.0 / sqrt(nx*nx + ny*ny);
    nx = nx * tmp;
    ny = ny * tmp;

    compute_T(q[i], nx, ny, T);

    r0 = T[ 0]*rhs[i][0] + T[ 1]*rhs[i][1] + T[ 2]*rhs[i][2] + T[ 3]*rhs[i][3];
    r1 = T[ 4]*rhs[i][0] + T[ 5]*rhs[i][1] + T[ 6]*rhs[i][2] + T[ 7]*rhs[i][3];
    r2 = T[ 8]*rhs[i][0] + T[ 9]*rhs[i][1] + T[10]*rhs[i][2] + T[11]*rhs[i][3];
    r3 = T[12]*rhs[i][0] + T[13]*rhs[i][1] + T[14]*rhs[i][2] + T[15]*rhs[i][3];

    rhs[i][0] = r0;
    rhs[i][1] = r1;
    rhs[i][2] = r2;
    rhs[i][3] = r3;
  }
  }
  
}

void invert_xi_eta(double (*rhs)[4], double (*q)[4], Dim *dim, Grid *grid){

  int i, j, k;
  double kx,ky,lx,ly, tmp;
  double N[16];
  double r0, r1, r2, r3;

  for(k=dim->nghost; k<dim->kmax + dim->nghost; k++){
  for(j=dim->nghost; j<dim->jmax + dim->nghost; j++){

    i = j*dim->jstride + k*dim->kstride;

    // "k" direction is k
    kx = grid->Sk[i][0];
    ky = grid->Sk[i][1];
    tmp = 1.0 / sqrt(kx*kx + ky*ky);
    kx = kx * tmp;
    ky = ky * tmp;

    // "l" direction is j
    lx = grid->Sj[i][0];
    ly = grid->Sj[i][1];
    tmp = 1.0 / sqrt(lx*lx + ly*ly);
    lx = lx * tmp;
    ly = ly * tmp;

    Tk_inv_Tl(kx,ky,lx,ly,N);

    r0 = N[ 0]*rhs[i][0] + N[ 1]*rhs[i][1] + N[ 2]*rhs[i][2] + N[ 3]*rhs[i][3];
    r1 = N[ 4]*rhs[i][0] + N[ 5]*rhs[i][1] + N[ 6]*rhs[i][2] + N[ 7]*rhs[i][3];
    r2 = N[ 8]*rhs[i][0] + N[ 9]*rhs[i][1] + N[10]*rhs[i][2] + N[11]*rhs[i][3];
    r3 = N[12]*rhs[i][0] + N[13]*rhs[i][1] + N[14]*rhs[i][2] + N[15]*rhs[i][3];

    rhs[i][0] = r0;
    rhs[i][1] = r1;
    rhs[i][2] = r2;
    rhs[i][3] = r3;
  }
  }
  
}

void xi_tridiagonal(double (*rhs)[4], double (*q)[4], double *dt,
		    double (*L)[4], double (*D)[4], double (*U)[4],
		    Dim *dim, Grid *grid){

  double kx, ky, jac, mag;
  double rho, irho, u, v, w, p, cn;
  double lambda_plus, lambda_minus, abs_lambda, lam0, lam2, lam3;
  double eps;
  int i, j, k;

  for(k=dim->nghost; k<dim->kmax + dim->nghost; k++){
  for(j=dim->nghost; j<dim->jmax + dim->nghost; j++){
  // for(k=1; k<dim->ktot-1; k++){
  // for(j=1; j<dim->jtot-1; j++){

    i   = j*dim->jstride + k*dim->kstride;

    jac = 1.0 / grid->V[i];

    // kx  = grid->Sj[i][0]*jac;
    // ky  = grid->Sj[i][1]*jac;
    kx  = 0.5*(grid->Sj[i][0] + grid->Sj[i+dim->jstride][0])*jac;
    ky  = 0.5*(grid->Sj[i][1] + grid->Sj[i+dim->jstride][1])*jac;

    mag = sqrt(kx*kx + ky*ky);

    rho    = q[i][0];
    irho   = 1.0 / rho;
    u      = q[i][1]*irho;
    v      = q[i][2]*irho;
    p      = (GAMMA-1.0)*(q[i][3] - 0.5*rho*(u*u + v*v));
    cn     = sqrt(GAMMA*p*irho)*mag;
    
    lam0 = kx*u + ky*v;
    lam2 = lam0 + cn;
    lam3 = lam0 - cn;

    eps = 0.08*mag;

    // lambda 0
    abs_lambda   = (1.05)*sqrt(lam0*lam0 + eps*eps);
    lambda_plus  = 0.5*(lam0 + abs_lambda);
    lambda_minus = 0.5*(lam0 - abs_lambda);

    L[i][0] =       (-lambda_plus               )*dt[i];
    D[i][0] = 1.0 + ( lambda_plus - lambda_minus)*dt[i];
    U[i][0] =       ( lambda_minus              )*dt[i];

    // lambda 1 = lambda 0
    L[i][1] =       (-lambda_plus               )*dt[i];
    D[i][1] = 1.0 + ( lambda_plus - lambda_minus)*dt[i];
    U[i][1] =       ( lambda_minus              )*dt[i];

    // lambda 2
    abs_lambda   = (1.05)*sqrt(lam2*lam2 + eps*eps);
    lambda_plus  = 0.5*(lam2 + abs_lambda);
    lambda_minus = 0.5*(lam2 - abs_lambda);

    L[i][2] =       (-lambda_plus               )*dt[i];
    D[i][2] = 1.0 + ( lambda_plus - lambda_minus)*dt[i];
    U[i][2] =       ( lambda_minus              )*dt[i];

    // lambda 3
    abs_lambda   = (1.05)*sqrt(lam3*lam3 + eps*eps);
    lambda_plus  = 0.5*(lam3 + abs_lambda);
    lambda_minus = 0.5*(lam3 - abs_lambda);

    L[i][3] =       (-lambda_plus               )*dt[i];
    D[i][3] = 1.0 + ( lambda_plus - lambda_minus)*dt[i];
    U[i][3] =       ( lambda_minus              )*dt[i];

  }
  }


  
  for(k=dim->nghost; k<dim->kmax + dim->nghost; k++){

    int start = dim->nghost*dim->jstride*4 + k*dim->kstride*4;

    solve_tridiagonal((double*)L, (double*)D, (double*)U, (double*)rhs,
		      dim->jmax, dim->jstride*4, start);
    solve_tridiagonal((double*)L, (double*)D, (double*)U, (double*)rhs,
		      dim->jmax, dim->jstride*4, start+1);
    solve_tridiagonal((double*)L, (double*)D, (double*)U, (double*)rhs,
		      dim->jmax, dim->jstride*4, start+2);
    solve_tridiagonal((double*)L, (double*)D, (double*)U, (double*)rhs,
		      dim->jmax, dim->jstride*4, start+3);    
  }

}

void eta_tridiagonal(double (*rhs)[4], double (*q)[4], double *dt,
		     double (*L)[4], double (*D)[4], double (*U)[4],
		     Dim *dim, Grid *grid){

  double kx, ky, jac, mag;
  double rho, irho, u, v, w, p, cn;
  double lambda_plus, lambda_minus, abs_lambda, lam0, lam2, lam3;
  double eps;
  int i, j, k;

  for(k=dim->nghost; k<dim->kmax + dim->nghost; k++){
  for(j=dim->nghost; j<dim->jmax + dim->nghost; j++){
  // for(j=1; j<dim->jtot-1; j++){
  // for(k=1; k<dim->ktot-1; k++){

    i   = j*dim->jstride + k*dim->kstride;

    jac = 1.0 / grid->V[i];

    kx  = 0.5*(grid->Sk[i][0] + grid->Sk[i+dim->kstride][0])*jac;
    ky  = 0.5*(grid->Sk[i][1] + grid->Sk[i+dim->kstride][1])*jac;

    mag = sqrt(kx*kx + ky*ky);

    rho    = q[i][0];
    irho   = 1.0 / rho;
    u      = q[i][1]*irho;
    v      = q[i][2]*irho;
    p      = (GAMMA-1.0)*(q[i][3] - 0.5*rho*(u*u + v*v));
    cn     = sqrt(GAMMA*p*irho)*mag;
    
    lam0 = kx*u + ky*v;
    lam2 = lam0 + cn;
    lam3 = lam0 - cn;

    eps = 0.08*mag;

    // lambda 0
    abs_lambda   = (1.05)*sqrt(lam0*lam0 + eps*eps);
    lambda_plus  = 0.5*(lam0 + abs_lambda);
    lambda_minus = 0.5*(lam0 - abs_lambda);

    L[i][0] =       (-lambda_plus               )*dt[i];
    D[i][0] = 1.0 + ( lambda_plus - lambda_minus)*dt[i];
    U[i][0] =       ( lambda_minus              )*dt[i];

    // lambda 1 = lambda 0
    L[i][1] =       (-lambda_plus               )*dt[i];
    D[i][1] = 1.0 + ( lambda_plus - lambda_minus)*dt[i];
    U[i][1] =       ( lambda_minus              )*dt[i];

    // lambda 2
    abs_lambda   = (1.05)*sqrt(lam2*lam2 + eps*eps);
    lambda_plus  = 0.5*(lam2 + abs_lambda);
    lambda_minus = 0.5*(lam2 - abs_lambda);

    L[i][2] =       (-lambda_plus               )*dt[i];
    D[i][2] = 1.0 + ( lambda_plus - lambda_minus)*dt[i];
    U[i][2] =       ( lambda_minus              )*dt[i];

    // lambda 3
    abs_lambda   = (1.05)*sqrt(lam3*lam3 + eps*eps);
    lambda_plus  = 0.5*(lam3 + abs_lambda);
    lambda_minus = 0.5*(lam3 - abs_lambda);

    L[i][3] =       (-lambda_plus               )*dt[i];
    D[i][3] = 1.0 + ( lambda_plus - lambda_minus)*dt[i];
    U[i][3] =       ( lambda_minus              )*dt[i];

  }
  }


  for(j=dim->nghost; j<dim->jmax + dim->nghost; j++){

    int start = dim->nghost*dim->kstride*4 + j*dim->jstride*4;

    solve_tridiagonal((double*)L, (double*)D, (double*)U, (double*)rhs,
		      dim->kmax, dim->kstride*4, start);
    solve_tridiagonal((double*)L, (double*)D, (double*)U, (double*)rhs,
		      dim->kmax, dim->kstride*4, start+1);
    solve_tridiagonal((double*)L, (double*)D, (double*)U, (double*)rhs,
		      dim->kmax, dim->kstride*4, start+2);
    solve_tridiagonal((double*)L, (double*)D, (double*)U, (double*)rhs,
		      dim->kmax, dim->kstride*4, start+3);
  }

}


void Euler::dadi(){

  double (*L)[4] = new double[dim->pts][4];
  double (*D)[4] = new double[dim->pts][4];
  double (*U)[4] = new double[dim->pts][4];
  
  // multiply rhs by T_xi_inv
  invert_xi(rhs, q, dim, grid);

  // Find L-D-U for Xi-direction, do tri-diagonal
  xi_tridiagonal(rhs, q, dt, L, D, U, dim, grid);

  // multiply rhs by T_eta_inv * T_xi
  invert_xi_eta(rhs, q, dim, grid);

  // Find L-D-U for Eta-direction, do tri-diagonal
  eta_tridiagonal(rhs, q, dt, L, D, U, dim, grid);
  
  // multiply rhs by T_eta
  invert_eta(rhs, q, dim, grid);

  delete L;
  delete D;
  delete U;
  
}
