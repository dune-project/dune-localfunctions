# set up project
project("dune-localfunctions" C CXX)

# general stuff
cmake_minimum_required(VERSION 3.13)

# guess dune-common build dir
if(NOT (dune-common_DIR OR dune-common_ROOT OR
        "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH}
  "${PROJECT_SOURCE_DIR}/cmake/modules")

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

add_subdirectory(doc)
add_subdirectory(dune)

if( DUNE_ENABLE_PYTHONBINDINGS )
  add_subdirectory(python)
  dune_python_install_package(PATH python)
endif()


finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
