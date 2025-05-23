# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
cmake_minimum_required(VERSION 3.16)

# set up project
project("dune-localfunctions" C CXX)

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

# deactivate global include-directories
dune_policy(SET DP_DEFAULT_INCLUDE_DIRS NEW)

# do not make all suggested dependencies required automatically
dune_policy(SET DP_SUGGESTED_MODULE_DEPENDENCIES_REQUIRED_DOWNSTREAM NEW)

# start a dune project with information from dune.module
dune_project()

# create library target and export it as Dune::LocalFunctions
dune_add_library(dunelocalfunctions INTERFACE
  EXPORT_NAME LocalFunctions
  LINK_LIBRARIES Dune::Geometry)

# set include directories to target
dune_default_include_directories(dunelocalfunctions INTERFACE)

add_subdirectory(doc)
add_subdirectory(dune)

if( DUNE_ENABLE_PYTHONBINDINGS )
  add_subdirectory(python)
endif()


finalize_dune_project()
