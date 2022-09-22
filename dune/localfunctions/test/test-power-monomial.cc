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
#include <dune/common/hybridutilities.hh>

#include <dune/localfunctions/meta/power.hh>
#include <dune/localfunctions/monomial.hh>

#include "geometries.hh"
#include "test-fe.hh"

// tolerance for floating-point comparisons
// be a little bit more lenient here than usual, the momom local basis
// is known to become more and more unstable with increasing order
static const double eps = 1e-8;
// stepsize for numerical differentiation
static const double delta = 1e-5;


template<int dimD, int dimR,int p>
static void Order(int &result)
{
        std::cout << "== Checking global-valued Power elements (with "
                  << "dimR=" << dimR << ") wrapping Monom elements (with "
                  << "dimD=" << dimD << ", p=" << p << ")" << std::endl;

        typedef TestGeometries<double, dimD> TestGeos;
        static const TestGeos testGeos;

        typedef typename TestGeos::Geometry Geometry;
        typedef Dune::MonomialFiniteElementFactory<Geometry, double, p>
        BackendFEFactory;
        BackendFEFactory backendFEFactory;

        typedef typename BackendFEFactory::FiniteElement BackendFE;
        Dune::PowerFiniteElementFactory<BackendFE, dimR> feFactory;

        for(std::size_t i = 0; i < testGeos.size(); ++i) {
          const Geometry &geo = testGeos[i];

          std::cout << "=== GeometryType " << geo.type() << std::endl;

          bool success =
            testFE(geo, feFactory.make(backendFEFactory.make(geo)), eps,
                   delta);

          if(success && result != 1)
            result = 0;
          else
            result = 1;
        }
}

template<int dimD,int dimR>
static void DimR(int &result)
{
   Dune::Hybrid::forEach(std::make_index_sequence<4>{},[&](auto i){Order<dimD,dimR,i>(result);});
}

template<int dimD>
static void DimD(int &result)
{
  Dune::Hybrid::forEach(std::make_index_sequence<4>{},[&](auto i){DimR<dimD,i>(result);});
}

int main(int argc, char** argv) {
  try {
    int result = 77;

    Dune::Hybrid::forEach(std::make_index_sequence<3>{},[&](auto i){DimD<i+1>(result);});

    return result;
  }
  catch (const Dune::Exception& e) {
    std::cerr << e << std::endl;
    throw;
  }
}
