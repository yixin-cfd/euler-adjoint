#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <stdio.h>
#include <boost/python.hpp>
//#define PY_ARRAY_UNIQUE_SYMBOL flow_ARRAY_API
//#include <numpy/ndarrayobject.h>
#include <string>
#include "grid.hpp"
#include "euler.hpp"
#include "meshgen.hpp"
#include "adadj.hpp"

BOOST_PYTHON_MODULE(libflow){

  using namespace boost::python;

  class_<Grid>("Grid", init<object,int>())
    .def("write_to_file", &Grid::write_to_file)
    ;

  class_<MeshGen>("MeshGen", init<object,int,double>())
    .def("write_to_file", &MeshGen::write_to_file)
    .def("poisson",&MeshGen::poisson)
    .def("get_mesh",&MeshGen::get_mesh)
    ;
  // class_<MeshGen>("MeshGen", init<object,int,double>());

  class_<Euler>("Euler", init<Grid*,std::string>())
    .def("say_hello", &Euler::say_hello)
    .def("write_solution", &Euler::write_solution)
    .def("go",&Euler::go)
    .def("take_steps",&Euler::take_steps)
    .def("pressure",&Euler::pressure)
    .def("Cl_Cd_Cm", &Euler::Cl_Cd_Cm)
    .def("save_restart", &Euler::save_restart)    
    .def("read_restart", &Euler::read_restart)
    ;

  class_<ADadj>("ADadj", init<Euler*>())
    .def("init", &ADadj::init)
    ;

}


