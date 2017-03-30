#ifdef __cplusplus
extern "C" {
#endif

void tmpflux_b(double *q_l, double *q_lb, double *q_r, double *q_rb, double *
	       rhs_m1, double *rhs_m1b, double *rhs, double *rhsb, Dim *dim, double S
	       [2], int mini, int maxi);

/* void tmpflux_b(double (*q)[4], double (*qb)[4], double (*rrhs)[4], double (*rrhsb)[4], */
/* 	       int lidx, int ridx, Dim *dim, double S[2], int mini, int maxi); */

/* void tmpflux_b(double *q_l, double *q_r, double *rhs_m1, double *rhs_m1b, */
/* 	       double *rhs, double *rhsb, Dim *dim, double xy1[2], double xy1b[2], */
/* 	       double xy2[2], double xy2b[2], int mini, int maxi); */



#ifdef __cplusplus
}
#endif
