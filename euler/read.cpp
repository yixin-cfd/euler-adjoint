#include "euler.hpp"
#include "yaml-cpp/yaml.h"


void Euler::read_inputs(std::string yaml_string){

  YAML::Node yaml = YAML::Load(yaml_string);

  inputs = new Inputs[1];
  
  int i;
  int nbc = yaml["bcs"].size();
  bc = new BC[nbc];

  inputs->M_inf = yaml["mach"].as<double>();
  inputs->aoa   = yaml["aoa"].as<double>();
  inputs->steps = yaml["steps"].as<int>();
  inputs->cfl   = yaml["cfl"].as<double>();
  inputs->ilhs  = yaml["ilhs"].as<int>();
  inputs->resid = yaml["resid"].as<int>();
  inputs->nbc   = nbc;

  for(i=0; i<nbc; i++){

    std::string face = yaml["bcs"][i]["face"].as<std::string>();
    std::string type = yaml["bcs"][i]["type"].as<std::string>();

    if      ( face.compare("jmin") == 0) bc[i].face = JMIN_FACE;
    else if ( face.compare("jmax") == 0) bc[i].face = JMAX_FACE;
    else if ( face.compare("kmin") == 0) bc[i].face = KMIN_FACE;
    else if ( face.compare("kmax") == 0) bc[i].face = KMAX_FACE;

    if      ( type.compare("periodic") == 0) bc[i].type = PERIODIC_BC;
    else if ( type.compare("farfield") == 0) bc[i].type = FARFIELD_BC;
    else if ( type.compare("wall"    ) == 0) bc[i].type = WALL_BC;

    bc[i].js = yaml["bcs"][i]["j"][0].as<int>();
    bc[i].je = yaml["bcs"][i]["j"][1].as<int>();
    bc[i].ks = yaml["bcs"][i]["k"][0].as<int>();
    bc[i].ke = yaml["bcs"][i]["k"][1].as<int>();

    //
    // Boundaries dont know about ghosts
    //
    if(bc[i].js > 0) bc[i].js += dim->nghost;
    if(bc[i].je > 0) bc[i].je += dim->nghost;
    if(bc[i].ks > 0) bc[i].ks += dim->nghost;
    if(bc[i].ke > 0) bc[i].ke += dim->nghost;

    //
    // Some boundaries are relative to the end dimension
    //
    bc[i].js = (bc[i].js + dim->jtot) % dim->jtot;
    bc[i].je = (bc[i].je + dim->jtot) % dim->jtot;
    bc[i].ks = (bc[i].ks + dim->ktot) % dim->ktot;
    bc[i].ke = (bc[i].ke + dim->ktot) % dim->ktot;

    switch(bc[i].face){
    case JMIN_FACE:
      bc[i].js = 0;
      bc[i].je = dim->nghost-1;
      break;
    case JMAX_FACE:
      bc[i].js = dim->jtot-dim->nghost;
      bc[i].je = dim->jtot-1;
      break;
    case KMIN_FACE:
      bc[i].ks = 0;
      bc[i].ke = dim->nghost-1;
      break;
    case KMAX_FACE:
      bc[i].ks = dim->ktot-dim->nghost;
      bc[i].ke = dim->ktot-1;
      break;
    }

  }

  // for(i=0; i<nbc; i++){
  //   printf("%d, %d: %d %d %d %d\n", (int)bc[i].face, (int)bc[i].type, 
  // 	   bc[i].js, bc[i].je, bc[i].ks, bc[i].ke);
  // }

}
