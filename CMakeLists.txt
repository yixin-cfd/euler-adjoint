#
# 2DFLOW CMakeLists.txt file
#
cmake_minimum_required(VERSION 2.8)

#
# Reset some options BEFORE declaring project
#
set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Project Build Type")

project(2DFLOW)
mark_as_advanced(FORCE CMAKE_INSTALL_PREFIX)

option(2DFLOW_DEBUG "Debugging Mode and Flags" OFF)


#
# Find Libraries
#
# YAML:
find_library(YAML_LIB yaml-cpp PATHS ${2DFLOW_SOURCE_DIR}/yaml/build NO_DEFAULT_PATH)
if(NOT YAML_LIB)
  message(FATAL_ERROR "YAML library not found, exiting.")
endif()

# Python & Boost
find_package(Boost COMPONENTS python REQUIRED)
find_package(PythonLibs REQUIRED)
message(STATUS "python libs: ${PYTHON_LIBRARIES}")
execute_process(COMMAND python -c "import mpi4py; print( mpi4py.get_include() )"
  OUTPUT_VARIABLE mpi4py_inc OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "mpi4py: ${mpi4py_inc}")
message(STATUS "boosts: ${Boost_LIBRARIES}")

#
# Include Directories
#
set(INCLUDE_DIRS ${2DFLOW_SOURCE_DIR}/include 
                 ${2DFLOW_SOURCE_DIR}/yaml/include
		 ${Boost_INCLUDE_DIRS}
		 ${PYTHON_INCLUDE_DIRS}
		 ${mpi4py_inc})

#
# Library: grid
#
file(GLOB grid_sources "${2DFLOW_SOURCE_DIR}/grid/*.cpp")

add_library(grid SHARED ${grid_sources})

set_target_properties(grid PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${2DFLOW_BINARY_DIR}/lib
  LIBRARY_OUTPUT_DIRECTORY ${2DFLOW_BINARY_DIR}/lib
  )
target_compile_definitions(grid PRIVATE ${DEBUG_DEFINITION})
target_include_directories(grid PUBLIC ${INCLUDE_DIRS})
target_link_libraries(grid ${YAML_LIB} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

#
# Library: euler
#
file(GLOB euler_sources "${2DFLOW_SOURCE_DIR}/euler/*.cpp")

add_library(euler SHARED ${euler_sources})

set_target_properties(euler PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${2DFLOW_BINARY_DIR}/lib
  LIBRARY_OUTPUT_DIRECTORY ${2DFLOW_BINARY_DIR}/lib
  )
target_compile_definitions(euler PRIVATE ${DEBUG_DEFINITION})
target_include_directories(euler PUBLIC ${INCLUDE_DIRS})
target_link_libraries(euler ${YAML_LIB} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

#
# Python Module for all libs
#
file(GLOB flow_sources   "${2DFLOW_SOURCE_DIR}/python_wrap/*.cpp" )

message(STATUS ${flow_sources})

add_library(flow SHARED ${flow_sources})

set_target_properties(flow PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${2DFLOW_BINARY_DIR}/lib
  LIBRARY_OUTPUT_DIRECTORY ${2DFLOW_BINARY_DIR}/lib
  )
target_compile_definitions(flow PRIVATE ${DEBUG_DEFINITION})
target_include_directories(flow PUBLIC ${INCLUDE_DIRS})
target_link_libraries(flow euler grid ${YAML_LIB} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
#target_link_libraries(flow ${YAML_LIB} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

