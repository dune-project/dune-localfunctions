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

#include "../lagrange/p0.hh"
#include "../lagrange/p1.hh"
#include "../lagrange/prismp1.hh"
#include "../lagrange/prismp2.hh"
#include "../lagrange/pyramidp1.hh"
#include "../lagrange/pyramidp2.hh"
#include "../lagrange/q1.hh"
#include "../lagrange/p23d.hh"
#include "../lagrange/pq22d.hh"
#include "../lagrange/pk.hh"
#include "../lagrange/qk.hh"

#include "../brezzidouglasmarini/brezzidouglasmarini1cube2d.hh"
#include "../brezzidouglasmarini/brezzidouglasmarini1cube3d.hh"
#include "../brezzidouglasmarini/brezzidouglasmarini2cube2d.hh"
#include "../brezzidouglasmarini/brezzidouglasmarini1simplex2d.hh"
#include "../brezzidouglasmarini/brezzidouglasmarini2simplex2d.hh"
#include "../dualmortarbasis.hh"
#include "../refined/refinedp1.hh"
#include "../refined/refinedp0.hh"
#include "../hierarchical/hierarchicalp2.hh"
#include "../hierarchical/hierarchicalp2withelementbubble.hh"
#include "../hierarchical/hierarchicalprismp2.hh"
#include "../rannacherturek/rannacherturek.hh"
#include "../raviartthomas/raviartthomassimplex.hh"
#include "../raviartthomas/raviartthomascube.hh"
#include "../monomial.hh"

#include "../common/virtualinterface.hh"
#include "../common/virtualwrappers.hh"
#include "../common/concepts.hh"
#include "../common/staticLoops.hh"
#include "../common/staticIf.hh"

#define FE_PAIR(...) std::make_pair( type<Dune:: __VA_ARGS__ >(), #__VA_ARGS__ )

template <class T>
struct Type { using type = T; };

template <class T>
Type<T> type() { return {}; }

template <class... Ts>
constexpr std::size_t size(std::tuple<Ts...> const&)
{
  return sizeof...(Ts);
}

// some dummy classes for false tests
struct A {};

struct B
{
  using Traits = double;
};

int main() try
{
  auto defaultTests = std::make_tuple(
    FE_PAIR( P0LocalFiniteElement<double,double,2> ),
    FE_PAIR( P1LocalFiniteElement<double,double,1> ),
    FE_PAIR( P1LocalFiniteElement<double,double,2> ),
    FE_PAIR( P1LocalFiniteElement<double,double,3> ),
    FE_PAIR( Q1LocalFiniteElement<double,double,1> ),
    FE_PAIR( Q1LocalFiniteElement<double,double,2> ),
    FE_PAIR( Q1LocalFiniteElement<double,double,3> ),
    FE_PAIR( PQ22DLocalFiniteElement<double,double> ),
    FE_PAIR( RefinedP1LocalFiniteElement<double,double,1> ),
    FE_PAIR( RefinedP1LocalFiniteElement<double,double,2> ),
    FE_PAIR( RefinedP1LocalFiniteElement<double,double,3> ),
    FE_PAIR( RefinedP0LocalFiniteElement<double,double,1> ),
    FE_PAIR( RefinedP0LocalFiniteElement<double,double,2> ),
    FE_PAIR( P23DLocalFiniteElement<double,double> ),
    FE_PAIR( HierarchicalP2LocalFiniteElement<double,double,1> ),
    FE_PAIR( HierarchicalP2LocalFiniteElement<double,double,2> ),
    FE_PAIR( HierarchicalP2LocalFiniteElement<double,double,3> ),
    FE_PAIR( HierarchicalPrismP2LocalFiniteElement<double,double> ),
    FE_PAIR( HierarchicalP2WithElementBubbleLocalFiniteElement<double,double,2> ),
    FE_PAIR( PrismP1LocalFiniteElement<double,double> ),
    FE_PAIR( PrismP2LocalFiniteElement<double,double> ),
    FE_PAIR( QkLocalFiniteElement<double,double,1,1> ),
    FE_PAIR( QkLocalFiniteElement<double,double,2,0> ),
    FE_PAIR( QkLocalFiniteElement<double,double,2,1> ),
    FE_PAIR( QkLocalFiniteElement<double,double,2,2> ),
    FE_PAIR( QkLocalFiniteElement<double,double,2,3> ),
    FE_PAIR( QkLocalFiniteElement<double,double,3,0> ),
    FE_PAIR( QkLocalFiniteElement<double,double,3,1> ),
    FE_PAIR( QkLocalFiniteElement<double,double,3,2> ),
    FE_PAIR( QkLocalFiniteElement<double,double,3,3> ),
    FE_PAIR( DualP1LocalFiniteElement<double,double,1> ),
    FE_PAIR( DualP1LocalFiniteElement<double,double,2> ),
    FE_PAIR( DualP1LocalFiniteElement<double,double,3> ),
    FE_PAIR( DualQ1LocalFiniteElement<double,double,1> ),
    FE_PAIR( DualQ1LocalFiniteElement<double,double,2> ),
    FE_PAIR( DualQ1LocalFiniteElement<double,double,3> ),
    FE_PAIR( RannacherTurekLocalFiniteElement<double,double,2> ),
    FE_PAIR( RannacherTurekLocalFiniteElement<double,double,3> ),
    FE_PAIR( BDM1Cube2DLocalFiniteElement<double,double> ),
    FE_PAIR( BDM2Cube2DLocalFiniteElement<double,double> ),
    FE_PAIR( BDM1Simplex2DLocalFiniteElement<double,double> ),
    FE_PAIR( RaviartThomasSimplexLocalFiniteElement<2,double,double> ),
    FE_PAIR( RaviartThomasSimplexLocalFiniteElement<2,double,double> ),
    FE_PAIR( RaviartThomasCubeLocalFiniteElement<double,double,2,0> ),
    FE_PAIR( RaviartThomasCubeLocalFiniteElement<double,double,3,0> ),
    FE_PAIR( RaviartThomasCubeLocalFiniteElement<double,double,2,1> ),
    FE_PAIR( RaviartThomasCubeLocalFiniteElement<double,double,3,1> ),
    FE_PAIR( RaviartThomasCubeLocalFiniteElement<double,double,2,2> ),
    FE_PAIR( BDM1Cube3DLocalFiniteElement<double,double> ),
    FE_PAIR( BDM2Simplex2DLocalFiniteElement<double,double> ),
    FE_PAIR( PyramidP1LocalFiniteElement<double,double> ),
    FE_PAIR( PyramidP2LocalFiniteElement<double,double> ),
    FE_PAIR( RaviartThomasCubeLocalFiniteElement<double,double,2,3> ),
    FE_PAIR( RaviartThomasCubeLocalFiniteElement<double,double,2,4> ),
    FE_PAIR( MonomialLocalFiniteElement<double,double,1,1> ),
    FE_PAIR( MonomialLocalFiniteElement<double,double,2,1> ),
    FE_PAIR( MonomialLocalFiniteElement<double,double,3,1> ),
    FE_PAIR( PkLocalFiniteElement<double,double,2,2> )
    );

  // run tests...
  bool success = Dune::staticAccumulate<0, size(defaultTests)>([&tests = defaultTests](auto const i)
  {
    using FE = typename std::decay_t< decltype( std::get<i>(tests).first ) >::type;

    bool b = Dune::Concept::isLocalFiniteElement<FE>();

    if (!b)
      std::cout << "concept(LocalFiniteElement) not satisfied for "
                << std::get<i>(tests).second << "\n";

    return b;
  }, true, [](bool a, bool b) { return a && b; });

  // ---------------------------------------------------------------------------
  // test virtualized FEs

  // run tests...
  success = Dune::staticAccumulate<0, size(defaultTests)>([&tests = defaultTests](auto const i)
  {
    using FE = typename std::decay_t< decltype( std::get<i>(tests).first ) >::type;

    using FE1 = Dune::LocalFiniteElementVirtualImp<FE>;
    bool b1 = Dune::Concept::isLocalFiniteElement<FE1>();

    using FE2 = Dune::LocalFiniteElementVirtualImp<Dune::LocalFiniteElementVirtualImp<FE>>;
    bool b2 = Dune::Concept::isLocalFiniteElement<FE2>();

    using FE3 = Dune::LocalFiniteElementVirtualInterface<typename FE::Traits::LocalBasisType::Traits>;
    bool b3 = Dune::Concept::isLocalFiniteElement<FE3>();

    if (!b1 || !b2 || !b3)
      std::cout << "concept(LocalFiniteElement) not satisfied for Virtual<"
                << std::get<i>(tests).second << ">\n";

    return b1 && b2 && b3;
  }, success, [](bool a, bool b) { return a && b; });


  // ---------------------------------------------------------------------------
  // false tests

  bool b0 = Dune::Concept::isLocalFiniteElement<A>();
  bool b1 = Dune::Concept::isLocalBasis<A>();
  bool b2 = Dune::Concept::isLocalCoefficients<A>();

  if (b0)
    std::cout << "concept(LocalFiniteElement) should not be satisfied for dummy class A\n";
  if (b1)
    std::cout << "concept(LocalBasis) should not be satisfied for dummy class A\n";
  if (b2)
    std::cout << "concept(LocalCoefficients) should not be satisfied for dummy class A\n";

  success = success && !(b0 || b1 || b2);

  b0 = Dune::Concept::isLocalFiniteElement<B>();
  b1 = Dune::Concept::isLocalBasis<B>();
  b2 = Dune::Concept::isLocalCoefficients<B>();

  if (b0)
    std::cout << "concept(LocalFiniteElement) should not be satisfied for dummy class B\n";
  if (b1)
    std::cout << "concept(LocalBasis) should not be satisfied for dummy class B\n";
  if (b2)
    std::cout << "concept(LocalCoefficients) should not be satisfied for dummy class B\n";

  success = success && !(b0 || b1 || b2);

  if (success)
    std::cout << "Success!\n";
  else
    std::cout << "Error!\n";

  return success ? 0 : 1;
}
catch (Dune::Exception e)
{
  std::cout << e << std::endl;
  return 1;
}
