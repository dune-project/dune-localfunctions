// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <utility>

#include <dune/common/hybridutilities.hh>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/lagrange/pqkfactory.hh>
#include <dune/localfunctions/dualmortarbasis/dualpq1factory.hh>
#include <dune/localfunctions/raviartthomas/raviartthomaslfecache.hh>

template<class FiniteElementCache>
static void test(Dune::GeometryType type)
{
  FiniteElementCache cache;

  using FiniteElement = typename FiniteElementCache::FiniteElementType;
  [[maybe_unused]] const FiniteElement& finiteElement = cache.get(type);
}

int main() {
  static constexpr std::size_t max_k = 3;
  Dune::Hybrid::forEach(std::make_index_sequence<max_k+1>{},[&](auto k)
          {
            constexpr int dim = 2;
            using FiniteElementCache = typename
                Dune::PQkLocalFiniteElementCache<double, double, dim, k>;
            test<FiniteElementCache>(Dune::GeometryTypes::simplex(dim));
            test<FiniteElementCache>(Dune::GeometryTypes::cube(dim));
          });
  {
    constexpr int dim = 2;
    using FiniteElementCache = typename
        Dune::DualPQ1LocalFiniteElementCache<double, double, dim>;
    test<FiniteElementCache>(Dune::GeometryTypes::simplex(dim));
    test<FiniteElementCache>(Dune::GeometryTypes::cube(dim));
  }

  {
    constexpr int dim = 2;
    constexpr int order = 0;
    using FiniteElementCache = typename
        Dune::RaviartThomasLocalFiniteElementCache<double, double, dim, order>;
    test<FiniteElementCache>(Dune::GeometryTypes::simplex(dim));
    test<FiniteElementCache>(Dune::GeometryTypes::cube(dim));
  }

  {
    constexpr int dim = 2;
    constexpr int order = 1;
    using FiniteElementCache = typename
        Dune::RaviartThomasLocalFiniteElementCache<double, double, dim, order>;
    test<FiniteElementCache>(Dune::GeometryTypes::simplex(dim));
    test<FiniteElementCache>(Dune::GeometryTypes::cube(dim));
  }

  {
    constexpr int dim = 3;
    constexpr int order = 0;
    using FiniteElementCache = typename
        Dune::RaviartThomasLocalFiniteElementCache<double, double, dim, order>;
    test<FiniteElementCache>(Dune::GeometryTypes::simplex(dim));
    test<FiniteElementCache>(Dune::GeometryTypes::cube(dim));
  }

  return 0;
}
