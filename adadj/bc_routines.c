#include "structures.h"

void periodic_bc(double (*q)[4], double (*rhs)[4], Dim *dim, BCface face, int j, int k){

  int idx, pidx;
  int pad = dim->nghost;

  idx = j*dim->jstride + k*dim->kstride;

  // find the periodic idx
  if(face == JMIN_FACE){
    pidx = (dim->jtot - 2*pad +j)*dim->jstride + k*dim->kstride;
  }
  if(face == JMAX_FACE){
    pidx = (j - dim->jtot + 2*pad)*dim->jstride + k*dim->kstride;
  }
  if(face == KMIN_FACE){
    pidx = j*dim->jstride + (dim->ktot - 2*pad + k)*dim->kstride;
  }
  if(face == KMAX_FACE){
    pidx = j*dim->jstride + (k - dim->ktot + 2*pad)*dim->kstride;
  }

  q[idx][0]   = q[pidx][0];
  q[idx][1]   = q[pidx][1];
  q[idx][2]   = q[pidx][2];
  q[idx][3]   = q[pidx][3];
  // rhs[idx][0]   = rhs[pidx][0];
  // rhs[idx][1]   = rhs[pidx][1];
  // rhs[idx][2]   = rhs[pidx][2];
  // rhs[idx][3]   = rhs[pidx][3];


}

void wall_bc(double (*q)[4], double (*rhs)[4], Dim *dim, 
	     double xy1[2], double xy2[2], int j, int k){

  int idx, widx, midx;
  int jstride = dim->jstride;
  int kstride = dim->kstride;
  int pad = dim->nghost;

  double Sx, Sy, imag2;

  Sx = - xy2[1] + xy1[1];
  Sy =   xy2[0] - xy1[0];

  double u, v, irho, v_dot_wall, p;

  idx = j*dim->jstride + k*dim->kstride;

  midx = j*jstride + (2*dim->nghost - k - 1)*kstride;
  widx = j*jstride + dim->nghost*kstride;

  // normalize wall vector
  imag2 = 1.0 / (Sx*Sx + Sy*Sy);

  q[idx][0] = q[midx][0];
  q[idx][1] = q[midx][1] - 2.0*Sx*Sx*imag2*q[midx][1] - 2.0*Sx*Sy*imag2*q[midx][2];
  q[idx][2] = q[midx][2] - 2.0*Sx*Sy*imag2*q[midx][1] - 2.0*Sy*Sy*imag2*q[midx][2];
  q[idx][3] = q[midx][3];


}

void farfield_bc(double q[4], Inputs *inputs){

  q[0] = inputs->rho_inf;
  q[1] = inputs->rho_inf * inputs->u_inf;
  q[2] = inputs->rho_inf * inputs->v_inf;
  q[3] = inputs->e_inf;

}
