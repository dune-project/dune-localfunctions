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

#include <dune/localfunctions/lagrange/q2.hh>

#include "geometries.hh"
#include "test-fe.hh"

template <int dim>
void test(const double& eps, const double& delta, int& result)
{
  std::cout << "== Checking global-valued Q2 " << dim << "D elements" << std::endl;

  Dune::GeometryType gt(Dune::GeometryTypes::cube(dim));

  typedef TestGeometries<double, dim> TestGeos;
  static const TestGeos testGeos;

  typedef typename TestGeos::Geometry Geometry;
  const Geometry &geo = testGeos.get(gt);

  Dune::Q2FiniteElementFactory<Geometry, double> feFactory;
  bool success = testFE(geo, feFactory.make(geo), eps, delta);

  if(success && result != 1)
    result = 0;
  else
    result = 1;
}

int main(int argc, char** argv) {
  try {
    // tolerance for floating-point comparisons
    static const double eps = 1e-9;
    // stepsize for numerical differentiation
    static const double delta = 1e-5;

    int result = 77;

    test<1>(eps, delta, result);
    test<2>(eps, delta, result);
    test<3>(eps, delta, result);

    return result;
  }
  catch (const Dune::Exception& e) {
    std::cerr << e << std::endl;
    throw;
  }
}
