#include "euler.hpp"

// #define DO_SIMPLE
#define EPS 0.125

void simpleflux(double* q_l, double* q_r, double* f, Dim *dim, double S[2]){

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


void roeflux(double* q_l, double* q_r, double* f, Dim *dim, double S[2]){

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

  double w1, w2, lambda_tilda;

  lambda_tilda = std::max(4.0*(eig_r1 - eig_l1),1e-12);
  w2 = 0.5*(1.0 + copysign(1.0,std::abs(eig_a1) - 0.5*lambda_tilda));
  w1 = 1.0 - w2;
  eig_a1 = w1*(eig_a1*eig_a1/lambda_tilda + 
  	       0.25*lambda_tilda ) + w2*( std::abs(eig_a1) );
  lambda_tilda = std::max(4.0*(eig_r2 - eig_r2),1e-12);
  w2 = 0.5*(1.0 + copysign(1.0,std::abs(eig_a2) - 0.5*lambda_tilda));
  w1 = 1.0 - w2;
  eig_a2 = w1*(eig_a2*eig_a2/lambda_tilda + 
  	       0.25*lambda_tilda ) + w2*( std::abs(eig_a2) );
  lambda_tilda = std::max(4.0*(eig_r3 - eig_l3),1e-12);
  w2 = 0.5*(1.0 + copysign(1.0,std::abs(eig_a3) - 0.5*lambda_tilda));
  w1 = 1.0 - w2;
  eig_a3 = w1*(eig_a3*eig_a3/lambda_tilda + 
  	       0.25*lambda_tilda ) + w2*( std::abs(eig_a3) );

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

}


void Euler::flux(){
  
  int j, k, idx, stride;

  //
  // J-direction
  //
  for(k=dim->nghost; k< dim->kmax+dim->nghost; k++){
  for(j=dim->nghost; j<=dim->jmax+dim->nghost; j++){

    idx = j*dim->jstride + k*dim->kstride;
#ifdef DO_SIMPLE
    simpleflux(q[idx-dim->jstride], q[idx], f[idx], dim, grid->Sj[idx]);
#else
    roeflux(q[idx-dim->jstride], q[idx], f[idx], dim, grid->Sj[idx]);
#endif
    

    rhs[idx-dim->jstride][0] += (f[idx-dim->jstride][0] - f[idx][0])*(j>dim->nghost);
    rhs[idx-dim->jstride][1] += (f[idx-dim->jstride][1] - f[idx][1])*(j>dim->nghost);
    rhs[idx-dim->jstride][2] += (f[idx-dim->jstride][2] - f[idx][2])*(j>dim->nghost);
    rhs[idx-dim->jstride][3] += (f[idx-dim->jstride][3] - f[idx][3])*(j>dim->nghost);

  }
  }

  //
  // K-direction
  //
  for(j=dim->nghost; j< dim->jmax+dim->nghost; j++){
  for(k=dim->nghost; k<=dim->kmax+dim->nghost; k++){

    idx = j*dim->jstride + k*dim->kstride;
#ifdef DO_SIMPLE    
    simpleflux(q[idx-dim->kstride], q[idx], f[idx], dim, grid->Sk[idx]);
#else    
    roeflux(q[idx-dim->kstride], q[idx], f[idx], dim, grid->Sk[idx]);
#endif    

    rhs[idx-dim->kstride][0] += (f[idx-dim->kstride][0] - f[idx][0])*(k>dim->nghost);
    rhs[idx-dim->kstride][1] += (f[idx-dim->kstride][1] - f[idx][1])*(k>dim->nghost);
    rhs[idx-dim->kstride][2] += (f[idx-dim->kstride][2] - f[idx][2])*(k>dim->nghost);
    rhs[idx-dim->kstride][3] += (f[idx-dim->kstride][3] - f[idx][3])*(k>dim->nghost);

  }
  }

}
