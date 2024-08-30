# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

add_subdirectory(p0)

install(FILES
  cache.hh
  emptypoints.hh
  equidistantpoints.hh
  interpolation.hh
  lagrangebasis.hh
  lagrangecoefficients.hh
  lagrangecube.hh
  lagrangelfecache.hh
  lagrangeprism.hh
  lagrangepyramid.hh
  lagrangesimplex.hh
  p0.hh
  pk1d.hh
  pk2d.hh
  pq22d.hh
  pqkfactory.hh
  q1.hh
  q2.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/lagrange)
