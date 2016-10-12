#include "grid.hpp"

void Grid::metrics(){

  int j, k, idx;
  int ktot = dim->ktot;
  int jtot = dim->jtot;
  int jstride = dim->jstride;
  int kstride = dim->kstride;

  //
  // "Face" Metrics
  //
  for(k=0; k<ktot-1; k++){
    for(j=0; j<jtot; j++){

      idx = j*jstride + k*kstride;

      Sj[idx][0] =   xy[idx+kstride][1] - xy[idx][1];
      Sj[idx][1] = - xy[idx+kstride][0] + xy[idx][0];

    }
  }

  for(j=0; j<jtot-1; j++){
    for(k=0; k<ktot; k++){

      idx = j*jstride + k*kstride;

      Sk[idx][0] = - xy[idx+jstride][1] + xy[idx][1];
      Sk[idx][1] =   xy[idx+jstride][0] - xy[idx][0];

    }
  }

  idx = (jtot-1)*jstride + (ktot-1)*kstride;

  //
  // Edge
  //
  k = ktot-1;
  for(j=0; j<jtot; j++){
    idx = j*jstride + k*kstride;

    Sj[idx][0] = Sj[idx-kstride][0];
    Sj[idx][1] = Sj[idx-kstride][1];
  }

  j = jtot-1;
  for(k=0; k<ktot; k++){
    idx = j*jstride + k*kstride;

    Sk[idx][0] = Sk[idx-jstride][0];
    Sk[idx][1] = Sk[idx-jstride][1];
  }

  //
  // Now Volume
  //
  int idx1, idx2, idx3, idx4;
  double diag1_x, diag1_y, diag2_x, diag2_y;
  for(k=dim->nghost; k<ktot-dim->nghost; k++){
    for(j=dim->nghost; j<jtot-dim->nghost; j++){
  // for(k=0; k<ktot-1; k++){
  //   for(j=0; j<jtot-1; j++){
      
      idx1 = j*jstride + k*kstride;
      idx2 = idx1 + jstride;
      idx3 = idx2 + kstride;
      idx4 = idx1 + kstride;

      diag1_x = xy[idx3][0]-xy[idx1][0];
      diag1_y = xy[idx3][1]-xy[idx1][1];
      diag2_x = xy[idx4][0]-xy[idx2][0];
      diag2_y = xy[idx4][1]-xy[idx2][1];

      V[idx1] = 0.5*(diag1_x*diag2_y - diag1_y*diag2_x);

    }
  }

  // //
  // // Edge
  // //
  // k = ktot-1;
  // for(j=0; j<jtot-1; j++){
  //   idx = j*jstride + k*kstride;
  //   V[idx] = V[idx-kstride];
  // }
  // j = jtot-1;
  // for(k=0; k<ktot; k++){ // plus last point too
  //   idx = j*jstride + k*kstride;
  //   V[idx] = V[idx-jstride];
  // }
  
  // 
  // Ghosts
  //
  // j-direction
  for(k=0; k < dim->ktot; k++){
    for(j=0; j<dim->nghost; j++){
      idx1 = dim->nghost*jstride + k*kstride;
      idx2 = j*jstride + k*kstride;
      V[idx2] = V[idx1];
    }
    for(j=dim->jtot-dim->nghost; j<dim->jtot; j++){
      idx1 = (dim->jtot-dim->nghost-1)*jstride + k*kstride;
      idx2 = j*jstride + k*kstride;
      V[idx2] = V[idx1];
    }    
  }
  // k-direction
  for(j=0; j<dim->jtot; j++){
    for(k=0; k<dim->nghost; k++){
      idx1 = j*jstride + dim->nghost*kstride;
      idx2 = j*jstride + k*kstride;
      V[idx2] = V[idx1];
    }
    for(k=dim->ktot-dim->nghost; k<dim->ktot; k++){
      idx1 = j*jstride + (dim->ktot-dim->nghost-1)*kstride;
      idx2 = j*jstride + k*kstride;
      V[idx2] = V[idx1];
    }
  }


}
