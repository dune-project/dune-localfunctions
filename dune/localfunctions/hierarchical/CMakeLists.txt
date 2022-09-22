# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

add_subdirectory(hierarchicalp2)
add_subdirectory(hierarchicalp2withelementbubble)
add_subdirectory(hierarchicalprismp2)

install(FILES
  hierarchicalp2.hh
  hierarchicalp2withelementbubble.hh
  hierarchicalprismp2.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/hierarchical)
