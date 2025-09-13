// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#include <dune/geometry/type.hh>

#include <dune/common/test/testsuite.hh>
#include <dune/localfunctions/lagrange/cache.hh>


int main ( int argc, char **argv )
{
  using namespace Dune;

  TestSuite tests;

  {
    DynamicLagrangeLocalFiniteElementCache<double,double,2> cache(3);
    const auto& fe1 = cache.get(Dune::GeometryTypes::simplex(2));
    const auto& fe2 = cache.get(Dune::GeometryTypes::cube(2));
    tests.check(fe1.localBasis().size() == 10);
    tests.check(fe2.localBasis().size() == 16);
  }
  {
    DynamicLagrangeLocalFiniteElementCache<double,double,3> cache(2);
    const auto& fe1 = cache.get(Dune::GeometryTypes::simplex(3));
    const auto& fe2 = cache.get(Dune::GeometryTypes::cube(3));
    const auto& fe3 = cache.get(Dune::GeometryTypes::prism);
    const auto& fe4 = cache.get(Dune::GeometryTypes::pyramid);
    tests.check(fe1.localBasis().size() == 10);
    tests.check(fe2.localBasis().size() == 27);
    tests.check(fe3.localBasis().size() == 18);
    tests.check(fe4.localBasis().size() == 14);
  }

  return tests.exit();
}
