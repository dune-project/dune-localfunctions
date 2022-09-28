# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

add_subdirectory(dualp1)
add_subdirectory(dualq1)
add_subdirectory(test)

install(FILES
  dualp1.hh
  dualq1.hh
  dualpq1factory.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/dualmortarbasis)
