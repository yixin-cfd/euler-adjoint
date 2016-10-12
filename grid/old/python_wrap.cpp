#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <stdio.h>
#include <boost/python.hpp>
#define PY_ARRAY_UNIQUE_SYMBOL grid_ARRAY_API
#include <numpy/ndarrayobject.h>
#include "grid.hpp"


BOOST_PYTHON_MODULE(libgrid){

  using namespace boost::python;

  import_array();

  class_<Grid>("Grid", init<object,int>())
    .def_readonly("jtot", &Grid::jtot)
    .def_readonly("ktot", &Grid::ktot)
    .def_readonly("jmax", &Grid::jmax)
    .def_readonly("kmax", &Grid::kmax)
    .def("write_to_file", &Grid::write_to_file)
    ;

}


