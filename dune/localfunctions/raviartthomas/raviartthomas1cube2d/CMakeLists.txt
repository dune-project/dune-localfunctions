# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

install(FILES
  raviartthomas1cube2dlocalbasis.hh
  raviartthomas1cube2dlocalcoefficients.hh
  raviartthomas1cube2dlocalinterpolation.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/raviartthomas/raviartthomas1cube2d)
