#include "euler.hpp"

template<BCface face>
void periodic_bc(BC bc, double (*q)[4], Dim *dim){

  int j, k, idx, pidx;
  int pad = dim->nghost;

  for(k=bc.ks; k<=bc.ke; k++){
  for(j=bc.js; j<=bc.je; j++){
      
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

      q[idx][0] = q[pidx][0];
      q[idx][1] = q[pidx][1];
      q[idx][2] = q[pidx][2];
      q[idx][3] = q[pidx][3];

  }
  }
}

template<BCface face>
void wall_bc(BC bc, double (*q)[4], Dim *dim, Grid* grid){

  int j, k, idx, widx, midx;
  int jstride = dim->jstride;
  int kstride = dim->kstride;
  int pad = dim->nghost;

  double Sx, Sy, imag, mag;

  double u, v, irho, v_dot_wall, p;

  for(k=bc.ks; k<=bc.ke; k++){
  for(j=bc.js; j<=bc.je; j++){
      
      idx = j*dim->jstride + k*dim->kstride;

      // find the periodic idx
      if (face==JMIN_FACE){
	midx = (2*dim->nghost - j - 1)*jstride + k*kstride;
	widx = dim->nghost*jstride + k*kstride;
	Sx = grid->Sj[widx][0];
	Sy = grid->Sj[widx][1];
      }
      if (face==KMIN_FACE){
	midx = j*jstride + (2*dim->nghost - k - 1)*kstride;
	widx = j*jstride + dim->nghost*kstride;
	Sx = grid->Sk[widx][0];
	Sy = grid->Sk[widx][1];
      }
      if (face==JMAX_FACE){
	midx = (2*(dim->jmax) + 2*dim->nghost - j - 1)*jstride + k*kstride;
	widx = (dim->nghost + (dim->jmax))*jstride + k*kstride;
	Sx = grid->Sj[widx][0];
	Sy = grid->Sj[widx][1];
      }
      if (face==KMAX_FACE){
	midx = j*jstride + (2*(dim->kmax) + 2*dim->nghost - k - 1)*kstride;
	widx = j*jstride + ((dim->kmax) + dim->nghost)*kstride;
	Sx = grid->Sk[widx][0];
	Sy = grid->Sk[widx][1];
      }

      // normalize wall vector
      mag = sqrt(Sx*Sx + Sy*Sy);
      imag = 1.0 / mag;
      // if(mag > 1.0e-30){
      // 	imag = 1.0 / mag;
      // } else {
      // 	imag = 0;
      // }
      
      Sx *= imag;
      Sy *= imag;

      irho = 1.0 / q[midx][0];

      // find velocity on other side of wall
      u = q[midx][1]*irho;
      v = q[midx][2]*irho;
      p = (GAMMA - 1.0)*(q[midx][3] - 0.5*q[midx][0]*(u*u + v*v));

      // reflect normal component
      v_dot_wall = Sx * u + Sy * v;
      u = u - 2.0*Sx*v_dot_wall;
      v = v - 2.0*Sy*v_dot_wall;

      q[idx][0] = q[midx][0];
      q[idx][1] = u*q[midx][0];
      q[idx][2] = v*q[midx][0];
      q[idx][3] = p / (GAMMA - 1.0) + 0.5*q[midx][0]*(u*u + v*v);

  }
  }

}

void farfield_bc(BC bc, double (*q)[4], Dim *dim, Inputs *inputs){

  int j, k, idx, pidx;
  int pad = dim->nghost;

  for(k=bc.ks; k<=bc.ke; k++){
  for(j=bc.js; j<=bc.je; j++){
      
      idx = j*dim->jstride + k*dim->kstride;

      q[idx][0] = inputs->rho_inf;
      q[idx][1] = inputs->rho_inf * inputs->u_inf;
      q[idx][2] = inputs->rho_inf * inputs->v_inf;
      q[idx][3] = inputs->e_inf;

  }
  }

}

// slow start
void ramp_bc(double ratio, BC bc, double (*q)[4], Dim *dim, Inputs *inputs){

  int j, k, idx, pidx;
  int pad = dim->nghost;

  for(k=bc.ks; k<=bc.ke; k++){
  for(j=bc.js; j<=bc.je; j++){
      
      idx = j*dim->jstride + k*dim->kstride;

      q[idx][0] = q[idx][0]*(ratio) + (1.0-ratio)*inputs->rho_inf;
      q[idx][1] = q[idx][1]*(ratio) + (1.0-ratio)*inputs->rho_inf * inputs->u_inf;
      q[idx][2] = q[idx][2]*(ratio) + (1.0-ratio)*inputs->rho_inf * inputs->v_inf;
      q[idx][3] = q[idx][3]*(ratio) + (1.0-ratio)*inputs->e_inf;

  }
  }

}


void Adjoint::boundary_conditions(){

  int ramp = (step_number < 30);
  double ratio = 1.0;

  for(int i=0; i<inputs->nbc; i++){

    switch(bc[i].face){
    case JMIN_FACE:
      if (bc[i].type == PERIODIC_BC) periodic_bc<JMIN_FACE>(bc[i], q, dim);
      if (bc[i].type == FARFIELD_BC) farfield_bc(bc[i], q, dim, inputs);
      if (bc[i].type == WALL_BC)     wall_bc<JMIN_FACE>(bc[i], q, dim, grid);
      break;
    case JMAX_FACE:
      if (bc[i].type == PERIODIC_BC) periodic_bc<JMAX_FACE>(bc[i], q, dim);
      if (bc[i].type == FARFIELD_BC) farfield_bc(bc[i], q, dim, inputs);
      if (bc[i].type == WALL_BC)     wall_bc<JMAX_FACE>(bc[i], q, dim, grid);
      break;
    case KMIN_FACE:
      if (bc[i].type == PERIODIC_BC) periodic_bc<KMIN_FACE>(bc[i], q, dim);
      if (bc[i].type == FARFIELD_BC) farfield_bc(bc[i], q, dim, inputs);
      if (bc[i].type == WALL_BC)     wall_bc<KMIN_FACE>(bc[i], q, dim, grid);
      break;
    case KMAX_FACE:
      if (bc[i].type == PERIODIC_BC) periodic_bc<KMAX_FACE>(bc[i], q, dim);
      if (bc[i].type == FARFIELD_BC) farfield_bc(bc[i], q, dim, inputs);
      if (bc[i].type == WALL_BC)     wall_bc<KMAX_FACE>(bc[i], q, dim, grid);
      break;

    }

    // // slow start
    // if(ramp && bc[i].type == WALL_BC){
    //   ratio = step_number * 1.0 / 30.0;
    //   ratio = (10.0 - 15.0*ratio + 6.0*ratio*ratio)*ratio*ratio*ratio;
    //   ramp_bc(ratio, bc[i], q, dim, inputs);
    // }
    
  }
}
