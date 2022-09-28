# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

install(FILES
  refinedp0localbasis.hh
  refinedp0localcoefficients.hh
  refinedp0localinterpolation.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/refined/refinedp0)
