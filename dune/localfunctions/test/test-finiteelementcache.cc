// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/common/hybridutilities.hh>
#include <dune/common/std/utility.hh>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/lagrange/pqkfactory.hh>
#include <dune/localfunctions/dualmortarbasis/dualpq1factory.hh>

template<class FiniteElementCache>
static void test(Dune::GeometryType type)
{
  FiniteElementCache cache;

  using FiniteElement = typename FiniteElementCache::FiniteElementType;
  DUNE_UNUSED const FiniteElement& finiteElement = cache.get(type);
}

int main() {
  static constexpr std::size_t max_k = 3;
  Dune::Hybrid::forEach(Dune::Std::make_index_sequence<max_k+1>{},[&](auto k)
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

  return 0;
}
