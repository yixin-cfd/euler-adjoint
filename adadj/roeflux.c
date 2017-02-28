#include "structures.h"

#define EPS 0.25

void jkroeflux(double* q_l, double* q_r, double* rhs_m1, double* rhs, Dim *dim,
	       double xy1[2], double xy2[2], int j_or_k, int is_j){

  double f[4];
  double S[2];

  S[0] =   xy2[1] - xy1[1];
  S[1] = - xy2[0] + xy1[0];

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

  double eig_l, eig_r;
  if(V_l > 0.0) 
    eig_l =  V_l + c_l;
  else
    eig_l = -V_l + c_l;
  if(V_r > 0.0) 
    eig_r =  V_r + c_r;
  else
    eig_r = -V_r + c_r;

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

  if(j_or_k > dim->nghost){
    rhs_m1[0] = rhs_m1[0] + f[0];
    rhs_m1[1] = rhs_m1[1] + f[1];
    rhs_m1[2] = rhs_m1[2] + f[2];
    rhs_m1[3] = rhs_m1[3] + f[3];
  }

  if(is_j==1 && j_or_k < dim->nghost + dim->jmax){
    rhs[0]    = rhs[0] - f[0];
    rhs[1]    = rhs[1] - f[1];
    rhs[2]    = rhs[2] - f[2];
    rhs[3]    = rhs[3] - f[3];
  }

  if(is_j==0 && j_or_k < dim->nghost + dim->kmax){
    rhs[0]    = rhs[0] - f[0];
    rhs[1]    = rhs[1] - f[1];
    rhs[2]    = rhs[2] - f[2];
    rhs[3]    = rhs[3] - f[3];
  }

}

void aajkroeflux(double* q_l, double* q_r, double* rhs_m1, double* rhs, Dim *dim,
	       double xy1[2], double xy2[2], int j_or_k, int is_j){

  double f[4];
  double S[2];

  S[0] =   xy2[1] - xy1[1];
  S[1] = - xy2[0] + xy1[0];

  //
  // Computing Roe Averages
  //
  double tmp;

  double rho_l  = q_l[0];
  double rho_r  = q_r[0];
  double rho_av = sqrt(rho_l*rho_r);
  double inv_rho_l  = 1.0/rho_l;
  double inv_rho_r  = 1.0/rho_r;

  double roe_wt_1 = sqrt(rho_l)/(sqrt(rho_l)+sqrt(rho_r));
  double roe_wt_2 = 1-roe_wt_1;

  double u_l  = q_l[1]*inv_rho_l;
  double u_r  = q_r[1]*inv_rho_r;
  double u_av = roe_wt_1*u_l + roe_wt_2*u_r;    

  double v_l  = q_l[2]*inv_rho_l;
  double v_r  = q_r[2]*inv_rho_r;
  double v_av = roe_wt_1*v_l + roe_wt_2*v_r;     

  double e_l  = q_l[3];
  double e_r  = q_r[3];

  double p_l  = (GAMMA-1.0)*(e_l -  0.5*rho_l*(u_l*u_l + v_l*v_l));
  double p_r  = (GAMMA-1.0)*(e_r -  0.5*rho_r*(u_r*u_r + v_r*v_r));

  double q_sqr_av = u_av*u_av + v_av*v_av;

  double h_l  = (e_l + p_l)*inv_rho_l;
  double h_r  = (e_r + p_r)*inv_rho_r;
  double h_av = roe_wt_1*h_l + roe_wt_2*h_r;     

  double c_l  = sqrt((GAMMA*p_l)*inv_rho_l);
  double c_r  = sqrt((GAMMA*p_r)*inv_rho_r);
  double c_av = sqrt((GAMMA-1.0)*(h_av - 0.5*(q_sqr_av))); 

  // face size and normals
  double mag = sqrt(S[0]*S[0] + S[1]*S[1]);
  double imag = 1.0/mag;
  double r1 = S[0]*imag;
  double r2 = S[1]*imag;

  // need vn grid here for lambdas below
  double V_l  = u_l *r1 + v_l *r2;
  double V_r  = u_r *r1 + v_r *r2;
  double V    = u_av*r1 + v_av*r2; 

  //
  // Computing Local Eigenvalues
  //
  double eig_a1 = V; 
  double eig_a2 = V + c_av; 
  double eig_a3 = V - c_av; 

  double eig_l1 = V_l; 
  double eig_l2 = V_l + c_l; 
  double eig_l3 = V_l - c_l; 

  double eig_r1 = V_r; 
  double eig_r2 = V_r + c_r; 
  double eig_r3 = V_r - c_r; 

  // Dylan: this looks like Harten's Entropy Correction but I'm not sure
  double lambda_tilda;
  eig_a1 = abs(eig_a1);
  eig_a2 = abs(eig_a2);
  eig_a3 = abs(eig_a3);
  //lambda_tilda = max(4.0*(eig_r1 - eig_l1) + 1e-6, 0.0);
  if(4.0*(eig_r1 - eig_l1) + 1e-6 > 0.0)
    lambda_tilda = 4.0*(eig_r1 - eig_l1) + 1e-6;
  else
    lambda_tilda = 0.0;
  if( eig_a1 < 0.5 * lambda_tilda ){
    eig_a1 = eig_a1*eig_a1 / (4.0*(eig_r1 - eig_l1) + 1e-6) +0.25*(4.0*(eig_r1 - eig_l1) + 1e-6);
  }

  //lambda_tilda = max(4.0*(eig_r2 - eig_l2) + 1e-6, 0.0);
  if(4.0*(eig_r2 - eig_l2) + 1e-6 > 0.0)
    lambda_tilda = 4.0*(eig_r1 - eig_l1) + 1e-6;
  else
    lambda_tilda = 0.0;
  if( eig_a2 < 0.5 * lambda_tilda ){
    eig_a2 = eig_a2*eig_a2 / (4.0*(eig_r2 - eig_l2) + 1e-6) +0.25*(4.0*(eig_r2 - eig_l2) + 1e-6);
  }

  //lambda_tilda = max(4.0*(eig_r3 - eig_l3) + 1e-6, 0.0);
  if(4.0*(eig_r3 - eig_l3) + 1e-6 > 0.0)
    lambda_tilda = 4.0*(eig_r1 - eig_l1) + 1e-6;
  else
    lambda_tilda = 0.0;

  if( eig_a3 < 0.5 * lambda_tilda ){
    eig_a3 = eig_a3*eig_a3 / (4.0*(eig_r3 - eig_l3) + 1e-6) +0.25*(4.0*(eig_r3 - eig_l3) + 1e-6);
  }

  double drho = rho_r - rho_l;
  double dp   = p_r - p_l;
  double du   = u_r - u_l;
  double dv   = v_r - v_l;
  double dV   = V_r - V_l;

  // quantity in parentheses blazek eq 4.90
  double drho_dp  = drho - dp/(c_av*c_av);
  // quantity in parentheses blazek eq 4.89, 4.91
  tmp = 1.0 / (2.0*c_av*c_av);
  double dp_p_rho_c = (dp + rho_av*c_av*dV)*tmp;
  double dp_m_rho_c = (dp - rho_av*c_av*dV)*tmp;

  double dF0, dF1, dF2, dF3; // delta F, blazek eqn. 4.89-4.91

  // eqn. 4.90
  dF0 = eig_a1 * (drho_dp);
  dF1 = eig_a1 * ( (drho_dp)*u_av         + rho_av*(du - dV*r1) );
  dF2 = eig_a1 * ( (drho_dp)*v_av         + rho_av*(dv - dV*r2) );
  dF3 = eig_a1 * ( (drho_dp)*q_sqr_av*0.5 + rho_av*(u_av*du + v_av*dv - V*dV ) );

  // eqn. 4.91
  dF0 = dF0 + eig_a2 * dp_p_rho_c;
  dF1 = dF1 + eig_a2 * dp_p_rho_c * (u_av + c_av*r1);
  dF2 = dF2 + eig_a2 * dp_p_rho_c * (v_av + c_av*r2);
  dF3 = dF3 + eig_a2 * dp_p_rho_c * (h_av + c_av*V);

  // eqn 4.89
  dF0 = dF0 + eig_a3 * dp_m_rho_c;
  dF1 = dF1 + eig_a3 * dp_m_rho_c * (u_av - c_av*r1);
  dF2 = dF2 + eig_a3 * dp_m_rho_c * (v_av - c_av*r2);
  dF3 = dF3 + eig_a3 * dp_m_rho_c * (h_av - c_av*V);

  double half_face = 0.5*mag;

  f[0] = half_face*( (rho_l*V_l             ) + (rho_r*V_r             ) - dF0);
  f[1] = half_face*( (rho_l*u_l*V_l + p_l*r1) + (rho_r*u_r*V_r + p_r*r1) - dF1);
  f[2] = half_face*( (rho_l*v_l*V_l + p_l*r2) + (rho_r*v_r*V_r + p_r*r2) - dF2);
  f[3] = half_face*( (e_l+p_l)*V_l + (e_r+p_r)*V_r                       - dF3);

  if(j_or_k > dim->nghost){
    rhs_m1[0] = rhs_m1[0] + f[0];
    rhs_m1[1] = rhs_m1[1] + f[1];
    rhs_m1[2] = rhs_m1[2] + f[2];
    rhs_m1[3] = rhs_m1[3] + f[3];
  }

  if(is_j==1 && j_or_k < dim->nghost + dim->jmax){
    rhs[0]    = rhs[0] - f[0];
    rhs[1]    = rhs[1] - f[1];
    rhs[2]    = rhs[2] - f[2];
    rhs[3]    = rhs[3] - f[3];
  }

  if(is_j==0 && j_or_k < dim->nghost + dim->kmax){
    rhs[0]    = rhs[0] - f[0];
    rhs[1]    = rhs[1] - f[1];
    rhs[2]    = rhs[2] - f[2];
    rhs[3]    = rhs[3] - f[3];
  }

  // rhs_m1[0] = rhs_m1[0] - f[0];
  // rhs_m1[1] = rhs_m1[1] - f[1];
  // rhs_m1[2] = rhs_m1[2] - f[2];
  // rhs_m1[3] = rhs_m1[3] - f[3];

  // rhs[0]    = rhs[0] + f[0];
  // rhs[1]    = rhs[1] + f[1];
  // rhs[2]    = rhs[2] + f[2];
  // rhs[3]    = rhs[3] + f[3];
  
}
