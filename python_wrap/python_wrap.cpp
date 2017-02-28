#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <stdio.h>
#include <boost/python.hpp>
//#define PY_ARRAY_UNIQUE_SYMBOL flow_ARRAY_API
//#include <numpy/ndarrayobject.h>
#include <string>
#include "grid.hpp"
#include "euler.hpp"
#include "slow_euler.hpp"
#include "meshgen.hpp"
#include "adadj.hpp"
#include "adjoint.hpp"

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

  class_<Slow_Euler>("Slow_Euler", init<Grid*,std::string>())
    .def("say_hello", &Slow_Euler::say_hello)
    .def("write_solution", &Slow_Euler::write_solution)
    .def("go",&Slow_Euler::go)
    .def("take_steps",&Slow_Euler::take_steps)
    .def("pressure",&Slow_Euler::pressure)
    .def("Cl_Cd_Cm", &Slow_Euler::Cl_Cd_Cm)
    .def("save_restart", &Slow_Euler::save_restart)    
    .def("read_restart", &Slow_Euler::read_restart)
    ;

  class_<ADadj>("ADadj", init<Euler*>())
    .def("init", &ADadj::init)
    .def("take_steps",&ADadj::take_steps)
    .def("save_restart", &ADadj::save_restart)    
    .def("read_restart", &ADadj::read_restart)
    .def("check",&ADadj::check)
    .def("sens_xd",&ADadj::sens_xd)
    ;

  class_<Adjoint>("Adjoint", init<Euler*>())
    .def("init", &Adjoint::init)
    .def("take_steps",&Adjoint::take_steps)
    .def("say_hello",&Adjoint::say_hello)
    ;

}


