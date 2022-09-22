# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

add_subdirectory(p0)

install(FILES
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
  p1.hh
  p23d.hh
  p2.hh
  pk1d.hh
  pk2d.hh
  pk3d.hh
  pk.hh
  pq22d.hh
  pqkfactory.hh
  prismp1.hh
  prismp2.hh
  pyramidp1.hh
  pyramidp2.hh
  q1.hh
  q2.hh
  qk.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/localfunctions/lagrange)
