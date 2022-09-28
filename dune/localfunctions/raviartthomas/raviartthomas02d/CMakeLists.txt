# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

install(FILES
  raviartthomas02dlocalbasis.hh
  raviartthomas02dlocalcoefficients.hh
  raviartthomas02dlocalinterpolation.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/raviartthomas/raviartthomas02d)
