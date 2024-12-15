# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

add_python_targets(localfunctions
  __init__
)

dune_add_pybind11_module(NAME _localfunctions)
target_link_libraries(_localfunctions PUBLIC Dune::Common Dune::Geometry Dune::LocalFunctions)
if(SKBUILD)
  install(TARGETS _localfunctions LIBRARY DESTINATION python/dune/localfunctions)
endif()
