#include "structures.h"

// J is cost
void pressure_cost(double q[4], Dim *dim, double *J, double p_desired){

  double p;

  // pressure at the wall is approximately the pressure at the cell right above the wall

  p = (GAMMA-1.0)*(q[3] - 0.5*(q[1]*q[1] + q[2]*q[2])/q[0]);

  J[0] = J[0] + abs(p - p_desired);

}

void lift_cost(double q[4], double xy1[2], double xy2[2], Dim *dim, double *J,
	       double upx, double upy){

  double p, Sx, Sy, Fy, Fx;
  Sx = - xy2[1] + xy1[1];
  Sy =   xy2[0] - xy1[0];
  
  // pressure at the wall is approximately the pressure at the cell right above the wall

  p = (GAMMA-1.0)*(q[3] - 0.5*(q[1]*q[1] + q[2]*q[2])/q[0]);

  Fy = -p*Sy;
  Fx = -p*Sx;

  J[0] = J[0] + Fx*upx + Fy*upy;


}
