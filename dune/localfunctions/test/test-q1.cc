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

#include <dune/localfunctions/lagrange/q1.hh>

#include "geometries.hh"
#include "test-fe.hh"

template<std::size_t dim>
void testQ1(int &result) {
  // tolerance for floating-point comparisons
  static const double eps = 1e-9;
  // stepsize for numerical differentiation
  static const double delta = 1e-5;

  std::cout << "== Checking global-valued Q1 elements (with dim=" << dim << ")"
            << std::endl;

  Dune::GeometryType gt(Dune::GeometryTypes::cube(dim));

  typedef TestGeometries<double, dim> TestGeos;
  static const TestGeos testGeos;

  typedef typename TestGeos::Geometry Geometry;
  const Geometry &geo = testGeos.get(gt);

  Dune::Q1FiniteElementFactory<Geometry, double> feFactory;
  bool success = testFE(geo, feFactory.make(geo), eps, delta);

  if(success && result != 1)
    result = 0;
  else
    result = 1;
}

int main(int argc, char** argv) {
  try {
    int result = 77;

    testQ1<0>(result);
    testQ1<1>(result);
    testQ1<2>(result);
    testQ1<3>(result);

    return result;
  }
  catch (const Dune::Exception& e) {
    std::cerr << e << std::endl;
    throw;
  }
}
