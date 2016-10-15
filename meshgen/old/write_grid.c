#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "cgnslib.h"

void init_cgns(){
  char basename[33];
  int icelldim,iphysdim,index_file,index_base;

  /* open CGNS file for write */
  if (cg_open("grid_c.cgns",CG_MODE_WRITE,&index_file)) cg_error_exit();

  /* create base (user can give any name) */
  sprintf(basename, "Base");
  icelldim=2;
  iphysdim=2;
  cg_base_write(index_file,basename,icelldim,iphysdim,&index_base);

  cg_close(index_file);

  // printf("Base index: %d\n", index_base);

}

void write_grid(int step, double *x, double *y, int jtot, int ktot){

  cgsize_t isize[3][2];
  int ni,nj,nk,i,j,k;
  int index_file;
  int index_zone,index_coord,index_base;
  char basename[33],zonename[33];
  //
  // open CGNS file for write
  //
  if (cg_open("grid_c.cgns",CG_MODE_MODIFY,&index_file)) cg_error_exit();

  // Base is already created, we want to create a zone

  // define zone name (user can give any name)
  sprintf(zonename, "Zone %02d", step);
  // vertex size
  isize[0][0]=jtot;
  isize[0][1]=ktot;
  // cell size
  isize[1][0]=isize[0][0]-1;
  isize[1][1]=isize[0][1]-1;
  // boundary vertex size (always zero for structured grids)
  isize[2][0]=0;
  isize[2][1]=0;

  index_base = 1; // always 1 for now

  //
  // Create zone
  //
  cg_zone_write(index_file,index_base,zonename,*isize,Structured,&index_zone);
  //
  // Write Grid Coordinates
  //
  cg_coord_write(index_file,index_base,index_zone,RealDouble,"CoordinateX",
		 x,&index_coord);
  cg_coord_write(index_file,index_base,index_zone,RealDouble,"CoordinateY",
		 y,&index_coord);
  /* close CGNS file */
  cg_close(index_file);

  // printf("\nSuccessfully wrote grid to file grid_c.cgns\n");

}
