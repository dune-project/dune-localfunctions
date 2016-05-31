// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstddef>
#include <iostream>
#include <cstdlib>
#include <vector>

#include <dune/common/function.hh>
#include <dune/localfunctions/common/staticLoops.hh>

#include "../lagrange/pyramidp1.hh"
#include "../lagrange/pyramidp2.hh"

#include "test-pyramid.hh"

template <class... Ts>
constexpr std::size_t size(std::tuple<Ts...> const&)
{
  return sizeof...(Ts);
}

int main() try
{
  using std::make_tuple;
  using std::make_pair;

  auto simplexGeometry = Dune::GeometryType(Dune::GeometryType::simplex, 2);

  // ---------------------------------------------------------------------------
  // Test some FiniteElements with disabled jacobian evaluation tests.
  // Test only whether evaluateGradient() matches with partial()

  auto disabledJacobianTests = make_tuple(
    /* 0*/ make_pair( Dune::PyramidP1LocalFiniteElement<double,double>(), "pyramidp1fem"),
    /* 1*/ make_pair( Dune::PyramidP2LocalFiniteElement<double,double>(), "pyramidp2fem")
    );

  // run tests...
  bool success = Dune::staticAccumulate<0, size(disabledJacobianTests)>([&tests = disabledJacobianTests](auto const i)
  {
    bool b = testFE(std::get<i>(tests).first);
    std::cout << "testFE(" << std::get<i>(tests).second << ") " << (b ? "succeeded\n" : "failed\n");
    return b;
  }, true, [](bool a, bool b) { return a && b; });

  return success ? 0 : 1;
}
catch (Dune::Exception e)
{
  std::cout << e << std::endl;
  return 1;
}
