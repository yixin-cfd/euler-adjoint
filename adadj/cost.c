#include "structures.h"

// J is cost
void pressure_cost(double q[4], Dim *dim, double *J, 
		   double p_desired){

  double p;

  // pressure at the wall is approximately the pressure at the cell right above the wall

  p = (GAMMA-1.0)*(q[3] - 0.5*(q[1]*q[1] + q[2]*q[2])/q[0]);

  J[0] = J[0] + abs(p - p_desired);

}
