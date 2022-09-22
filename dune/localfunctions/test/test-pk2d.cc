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
#include <utility>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/hybridutilities.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/generalvertexorder.hh>

#include <dune/localfunctions/lagrange/pk2d.hh>

#include "geometries.hh"
#include "test-fe.hh"

// tolerance for floating-point comparisons
static const double eps = 1e-9;
// stepsize for numerical differentiation
static const double delta = 1e-5;

template<int k>
static void test(int &result)
{
    std::cout << "== Checking global-valued Pk2D elements (with k=" << k << ")"
              << std::endl;

    Dune::GeometryType gt(Dune::GeometryTypes::triangle);

    typedef TestGeometries<double, 2> TestGeos;
    static const TestGeos testGeos;

    typedef TestGeos::Geometry Geometry;
    const Geometry &geo = testGeos.get(gt);

    std::size_t vertexIds[] = {0, 1, 2};
    Dune::GeneralVertexOrder<2, std::size_t>
    vo(gt, vertexIds+0, vertexIds+3);

    Dune::Pk2DFiniteElementFactory<Geometry, double, k> feFactory;
    bool success = testFE(geo, feFactory.make(geo, vo), eps, delta);

    if(success && result != 1)
      result = 0;
    else
      result = 1;
}

int main(int argc, char** argv) {
  try {
    int result = 77;

    static constexpr std::size_t max_k = 20;
    Dune::Hybrid::forEach(std::make_index_sequence<max_k+1>{},[&](auto i){test<i>(result);});

    return result;
  }
  catch (const Dune::Exception& e) {
    std::cerr << e << std::endl;
    throw;
  }
}
