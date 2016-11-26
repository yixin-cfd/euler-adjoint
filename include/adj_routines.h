#include "structures.h"

#ifdef __cplusplus
extern "C" {
#endif


  void pressure_cost(double q[4], Dim *dim, double *J, double p_desired);

  void pressure_cost_b(double q[4], double qb[4], Dim *dim, double *J, double *Jb, double p_desired);



#ifdef __cplusplus
}
#endif
