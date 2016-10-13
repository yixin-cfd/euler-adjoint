#ifndef STRUCTURES_H
#define STRUCTURES_H
#include "constants.h"

#define GAMMA 1.4

enum BCface { JMIN_FACE, JMAX_FACE, KMIN_FACE, KMAX_FACE };

enum BCtype { WALL_BC, FARFIELD_BC, PERIODIC_BC };

struct Inputs {
  double M_inf;
  int steps;
  double cfl;
  double aoa;
  double u_inf, v_inf, p_inf, e_inf, rho_inf;
  int nbc, resid, ilhs;
};

struct BC {
  int js, ks, je, ke;
  BCface face;
  BCtype type;
};

struct Dim {
  int jmax, kmax;
  int jtot, ktot, nghost;
  int pts;
  int jstride, kstride;
};
#endif
