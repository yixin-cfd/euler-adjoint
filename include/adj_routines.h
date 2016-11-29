#include "structures.h"

#ifdef __cplusplus
extern "C" {
#endif


  void pressure_cost(double q[4], Dim *dim, double *J, double p_desired);
  void timestep(double q[4], double Sj[2], double Sk[2],  
		double *V, double cfl, double *dt);

  void pressure_cost_b(double q[4], double qb[4], Dim *dim, double *J, double *Jb, double p_desired);
  void wall_bc_b(double (*q)[4], double (*qb)[4], Dim *dim, double xy1[2], 
	       double xy1b[2], double xy2[2], double xy2b[2], int j, int k);
  void periodic_bc_b(double (*q)[4], double (*qb)[4], Dim *dim, BCface face, int
		     j, int k);

#ifdef __cplusplus
}
#endif
