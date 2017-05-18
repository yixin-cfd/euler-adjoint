# Adjoint Euler Solver #

This is a year-long project for AMSC 663-664, Fall 2016-Spring 2017.

## Main Source Directories ##

*meshgen     - 2D Elliptic Mesh Generation Library
*grid        - basic grid utilities for a 2D CFD grid
*euler       - Euler flow solver
*adadj       - Auto-Differentiated Adjoint
*adjoint     - By-Hand Adjoint
*python_wrap - Boost-Python module to make C++ accessible to Python
*python      - Python helper functions and classes

## Additional Directories ##

*sample_run_directories
*slow_euler  - A debug version of the Euler code with some routines re-written for auto-differentiation
*include     - include directory with all headers
*yaml        - External C++ library to parse YAML input files

# Building and Running #

The Adjoint-Euler code uses cmake to build all source files. This is intended for an "out of source build", namely used like this:

```
#!bash
mkdir build
cd build
cmake ..
make
```

The major dependency is Boost-Python, which on linux is as simple as running:

```
#!bash
sudo apt-get install libboost-python-dev
```

On mac with the Homebrew package manager in **theory** ( in my version of Homebrew the boost-python package was missing some binaries. See [here](http://www.pyimagesearch.com/2015/04/27/installing-boost-and-boost-python-on-osx-with-homebrew/)):

```
#!bash
brew install boost --with-python
brew install boost-python
```

I am not at all familiar with Windows for Python or C++ development.