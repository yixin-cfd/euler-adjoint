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
message(STATUS "boosts: ${Boost_LIBRARIES}")

#
# Include Directories
#
set(INCLUDE_DIRS ${2DFLOW_SOURCE_DIR}/include 
                 ${2DFLOW_SOURCE_DIR}/yaml/include
		 ${Boost_INCLUDE_DIRS}
		 ${PYTHON_INCLUDE_DIRS}
		 ${mpi4py_inc})

file(GLOB python_utils "${2DFLOW_SOURCE_DIR}/python_wrap/python_hepers.cpp")

#
# Library: meshgen
#
file(GLOB meshgen_sources "${2DFLOW_SOURCE_DIR}/meshgen/*.cpp" "${2DFLOW_SOURCE_DIR}/meshgen/*.c")

add_library(meshgen SHARED ${meshgen_sources})

set_property( TARGET meshgen APPEND_STRING PROPERTY COMPILE_FLAGS -Wall )

set_target_properties(meshgen PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${2DFLOW_BINARY_DIR}/lib
  LIBRARY_OUTPUT_DIRECTORY ${2DFLOW_BINARY_DIR}/lib
  )
target_compile_definitions(meshgen PRIVATE ${DEBUG_DEFINITION})
target_include_directories(meshgen PUBLIC ${INCLUDE_DIRS})
target_link_libraries(meshgen ${YAML_LIB} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
	       
#
# Library: grid
#
file(GLOB grid_sources "${2DFLOW_SOURCE_DIR}/grid/*.cpp" "${2DFLOW_SOURCE_DIR}/python_wrap/python_hepers.cpp")

add_library(grid SHARED ${grid_sources})

set_property( TARGET grid APPEND_STRING PROPERTY COMPILE_FLAGS -Wall )

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
file(GLOB euler_sources "${2DFLOW_SOURCE_DIR}/euler/*.c" "${2DFLOW_SOURCE_DIR}/euler/*.cpp" "${2DFLOW_SOURCE_DIR}/python_wrap/python_hepers.cpp")

add_library(euler SHARED ${euler_sources})

set_property( TARGET euler APPEND_STRING PROPERTY COMPILE_FLAGS -Wall )

set_target_properties(euler PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${2DFLOW_BINARY_DIR}/lib
  LIBRARY_OUTPUT_DIRECTORY ${2DFLOW_BINARY_DIR}/lib
  )
target_compile_definitions(euler PRIVATE ${DEBUG_DEFINITION})
target_include_directories(euler PUBLIC ${INCLUDE_DIRS})
target_link_libraries(euler ${YAML_LIB} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

#
# Library: slow_euler
#
file(GLOB slow_euler_sources "${2DFLOW_SOURCE_DIR}/slow_euler/*.c" "${2DFLOW_SOURCE_DIR}/slow_euler/*.cpp" "${2DFLOW_SOURCE_DIR}/python_wrap/python_hepers.cpp")

add_library(slow_euler SHARED ${slow_euler_sources})

set_property( TARGET slow_euler APPEND_STRING PROPERTY COMPILE_FLAGS -Wall )

set_target_properties(slow_euler PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${2DFLOW_BINARY_DIR}/lib
  LIBRARY_OUTPUT_DIRECTORY ${2DFLOW_BINARY_DIR}/lib
  )
target_compile_definitions(slow_euler PRIVATE ${DEBUG_DEFINITION})
target_include_directories(slow_euler PUBLIC ${INCLUDE_DIRS})
target_link_libraries(slow_euler ${YAML_LIB} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})


#
# Library: ADadj
#
file(GLOB adadj_sources "${2DFLOW_SOURCE_DIR}/adadj/*.c" "${2DFLOW_SOURCE_DIR}/adadj/*.cpp" "${2DFLOW_SOURCE_DIR}/python_wrap/python_hepers.cpp")

add_library(adadj SHARED ${adadj_sources})

set_property( TARGET adadj APPEND_STRING PROPERTY COMPILE_FLAGS -Wall )

set_target_properties(adadj PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${2DFLOW_BINARY_DIR}/lib
  LIBRARY_OUTPUT_DIRECTORY ${2DFLOW_BINARY_DIR}/lib
  )

list(APPEND INCLUDE_DIRS "$ENV{TAPENADE_HOME}/include")

target_compile_definitions(adadj PRIVATE ${DEBUG_DEFINITION})
target_include_directories(adadj PUBLIC ${INCLUDE_DIRS})
target_link_libraries(adadj ${YAML_LIB} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

#
# Library: Adjoint
#
file(GLOB adjoint_sources "${2DFLOW_SOURCE_DIR}/adjoint/*.cpp" "${2DFLOW_SOURCE_DIR}/python_wrap/python_hepers.cpp")

add_library(adjoint SHARED ${adjoint_sources})

set_property( TARGET adjoint APPEND_STRING PROPERTY COMPILE_FLAGS -Wall )

set_target_properties(adjoint PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${2DFLOW_BINARY_DIR}/lib
  LIBRARY_OUTPUT_DIRECTORY ${2DFLOW_BINARY_DIR}/lib
  )

list(APPEND INCLUDE_DIRS "$ENV{TAPENADE_HOME}/include")

target_compile_definitions(adjoint PRIVATE ${DEBUG_DEFINITION})
target_include_directories(adjoint PUBLIC ${INCLUDE_DIRS})
target_link_libraries(adjoint ${YAML_LIB} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})


#
# Python Module for all libs
#
file(GLOB flow_sources "${2DFLOW_SOURCE_DIR}/python_wrap/*.cpp" )

message(STATUS ${flow_sources})

add_library(flow SHARED ${flow_sources})

set_target_properties(flow PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${2DFLOW_BINARY_DIR}/lib
  LIBRARY_OUTPUT_DIRECTORY ${2DFLOW_BINARY_DIR}/lib
  )
target_compile_definitions(flow PRIVATE ${DEBUG_DEFINITION})
target_include_directories(flow PUBLIC ${INCLUDE_DIRS})
target_link_libraries(flow euler grid meshgen adadj adjoint slow_euler ${YAML_LIB} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
#target_link_libraries(flow ${YAML_LIB} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

