# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

install(FILES
  raviartthomas12dlocalbasis.hh
  raviartthomas12dlocalcoefficients.hh
  raviartthomas12dlocalinterpolation.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/raviartthomas/raviartthomas12d)
