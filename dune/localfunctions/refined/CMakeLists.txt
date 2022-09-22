# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

add_subdirectory(common)
add_subdirectory(refinedp0)
add_subdirectory(refinedp1)

install(FILES
  refinedp0.hh
  refinedp1.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/refined)
