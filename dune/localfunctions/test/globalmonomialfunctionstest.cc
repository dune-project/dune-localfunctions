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
#include <dune/common/hybridutilities.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/generalvertexorder.hh>

#include <dune/localfunctions/monomial.hh>

#include "geometries.hh"
#include "test-fe.hh"

// tolerance for floating-point comparisons
const double eps = 1e-9;
// stepsize for numerical differentiation
const double delta = 1e-5;

template<int dim,int p>
static void Order(int &result)
{
      std::cout << "== Checking global-valued monomial elements (with "
                << "dim=" << dim << ", p=" << p << ")" << std::endl;

      typedef TestGeometries<double, dim> TestGeos;
      static const TestGeos testGeos;

      typedef typename TestGeos::Geometry Geometry;
      Dune::MonomialFiniteElementFactory<Geometry, double, p> feFactory;

      for(std::size_t i = 0; i < testGeos.size(); ++i) {
        const Geometry &geo = testGeos[i];

        std::cout << "=== GeometryType " << geo.type() << std::endl;

        bool success = testFE(geo, feFactory.make(geo), eps, delta);

        if(success && result != 1)
          result = 0;
        else
          result = 1;
      }
}

template<int dim>
static void Dim(int &result)
{
  Dune::Hybrid::forEach(std::make_index_sequence<4>{},[&](auto i){Order<dim,i>(result);});
}

int main(int argc, char** argv) {
  try {
    int result = 77;

    Dune::Hybrid::forEach(std::make_index_sequence<3>{},[&](auto i){Dim<i+1>(result);});

    return result;
  }
  catch (const Dune::Exception& e) {
    std::cerr << e << std::endl;
    throw;
  }
}
