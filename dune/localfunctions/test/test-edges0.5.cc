// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstddef>
#include <iostream>
#include <ostream>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/generalvertexorder.hh>

#include <dune/localfunctions/whitney/edges0.5.hh>

#include "geometries.hh"
#include "test-fe.hh"

template<std::size_t dim>
void testEdgeS0_5(int &result) {
  // tolerance for floating-point comparisons
  static const double eps = 1e-9;
  // stepsize for numerical differentiation
  static const double delta = 1e-5;

  std::cout << "== Checking global-valued EdgeS0_5 elements (with "
            << "dim=" << dim << ")" << std::endl;

  Dune::GeometryType gt(Dune::GeometryTypes::simplex(dim));

  typedef TestGeometries<double, dim> TestGeos;
  static const TestGeos testGeos;

  typedef typename TestGeos::Geometry Geometry;
  const Geometry &geo = testGeos.get(gt);

  static_assert(dim <= 3, "Need to update vertexIds array for dim > 3");
  std::size_t vertexIds[] = {0, 1, 2, 3};
  Dune::GeneralVertexOrder<dim, std::size_t>
  vo(gt, vertexIds+0, vertexIds+dim+1);

  Dune::EdgeS0_5FiniteElementFactory<Geometry, double> feFactory;
  bool success = testFE(geo, feFactory.make(geo, vo), eps, delta);

  if(success && result != 1)
    result = 0;
  else
    result = 1;
}

int main(int argc, char** argv) {
  try {
    int result = 77;

    testEdgeS0_5<2>(result);
    testEdgeS0_5<3>(result);

    return result;
  }
  catch (const Dune::Exception& e) {
    std::cout << e << std::endl;
    throw;
  }
}
