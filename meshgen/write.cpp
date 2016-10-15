#include "meshgen.hpp"

int MeshGen::write_to_file(std::string s){

  FILE *fid;
  int pts = dim->pts;

  try {
    fid = fopen(s.c_str(),"w");
  } catch(int error){
    printf("*** Could not open file for writing\n");
    return 1;
  }

  int count = 0;
  fprintf(fid, "%d %d\n", dim->jtot, dim->ktot);

  for(int i=0; i<pts; i++){
    fprintf(fid, "%25.16e ", x[i]);
    count++;
    if(count%2 == 0)
      fprintf(fid, "\n");
  }
  for(int i=0; i<pts; i++){
    fprintf(fid, "%25.16e ", y[i]);
    count++;
    if(count%2 == 0)
      fprintf(fid, "\n");
  }

  fclose(fid);

  return 0;

}
