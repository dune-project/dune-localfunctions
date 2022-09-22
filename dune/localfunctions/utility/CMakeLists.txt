# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

install(FILES
  basisevaluator.hh
  basismatrix.hh
  basisprint.hh
  coeffmatrix.hh
  defaultbasisfactory.hh
  dglocalcoefficients.hh
  field.hh
  interpolationhelper.hh
  l2interpolation.hh
  lfematrix.hh
  localfiniteelement.hh
  monomialbasis.hh
  multiindex.hh
  polynomialbasis.hh
  tensor.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/utility)
