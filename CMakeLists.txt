# set up project
project("dune-localfunctions" C CXX)

# general stuff
cmake_minimum_required(VERSION 2.8.6)

#find dune-common and set the module path
find_package(dune-common)
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH}
  "${CMAKE_SOURCE_DIR}/cmake/modules")

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

# $Id: duneproject 5232 2008-07-03 09:34:11Z christi $

# we need the module file to be able to build via dunecontrol
set(EXTRA_DIST dune.module)

# don't follow the full GNU-standard
# we need automake 1.9 or newer
set(AUTOMAKE_OPTIONS  foreign 1.9)

set(SUBDIRS  dune doc m4)

foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt "EXCLUDE_FROM_ALL")
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})

finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
