#include "structures.h"

#ifdef __cplusplus
extern "C" {
#endif

  void pressure_cost(double q[4], Dim *dim, double *J, double p_desired, Inputs* inputs);
  void lift_cost(double q[4], double xy1[2], double xy2[2], Dim *dim, double *J, double ux, double uy);
  void ad_timestep(double q[4], double xy1[2], double xy2[2],  double xy3[2],
		   double cfl, double *dt);

  void dadi_b(double *q, double *qb, double *rhs, double *rhsb, double *dt, 
	      double (*Sj)[2], double (*Sk)[2], double *V, Dim *dim);

  //
  // QDIFF routines
  //
  void pressure_cost_b(double q[4], double qb[4], Dim *dim, double *J, double *
		       Jb, double p_desired, Inputs* inputs);

  void ad_timestep_b(double q[4], double qb[4], double xy1[2], double xy2[2], 
		     double xy3[2], double cfl, double *dt, double *dtb);

  void periodic_bc_b(double (*q)[4], double (*qb)[4], double (*rhs)[4], 
		     Dim *dim, BCface face, int j, int k);

  void wall_bc_b(double (*q)[4], double (*qb)[4], double (*rhs)[4], 
		 Dim *dim, double xy1[2], double xy2[2], int j, int k);


  void jkroeflux_b(double *q_l, double *q_lb, double *q_r, double *q_rb, double *
		   rhs_m1, double *rhs_m1b, double *rhs, double *rhsb, Dim *dim, double 
		   xy1[2], double xy2[2], int j_or_k, int is_j);

  void lift_cost_b(double q[4], double qb[4], double xy1[2], double xy2[2], Dim 
		   *dim, double *J, double *Jb, double upx, double upy);


  //
  // XDIFF routines
  //
  void ad_timestep_bx(double q[4], double qb[4], double xy1[2], double xy1b[2], 
		   double xy2[2], double xy2b[2], double xy3[2], double xy3b[2], double 
		   cfl, double *dt, double *dtb);

  void wall_bc_bx(double (*q)[4], double (*qb)[4], double (*rhs)[4], 
		  Dim *dim, double xy1[2], double xy1b[2], double xy2[2], double xy2b[2], 
		  int j, int k);

  void jkroeflux_bx(double *q_l, double *q_lb, double *q_r, double *q_rb, double 
		    *rhs_m1, double *rhs_m1b, double *rhs, double *rhsb, Dim *dim, double 
		    xy1[2], double xy1b[2], double xy2[2], double xy2b[2], int j_or_k, int is_j);

  void lift_cost_bx(double q[4], double qb[4], double xy1[2], double xy1b[2], 
		    double xy2[2], double xy2b[2], Dim *dim, double *J, double *Jb,
		    double upx, double upy);



#ifdef __cplusplus
}
#endif
