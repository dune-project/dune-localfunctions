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
#include <dune/localfunctions/dualmortarbasis.hh>
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

#include "test-localfe.hh"

// tmp for testing arbitrary order finite elements
template<unsigned int d, int k>
struct PkLocalFiniteElementTest
{
  static bool test()
  {
    bool success = true;

    Dune::PkLocalFiniteElement<double,double,d,k> pklfem;
    TEST_FE(pklfem);

    return PkLocalFiniteElementTest<d, k-1>::test() and success;
  }
};

template<unsigned int d>
struct PkLocalFiniteElementTest<d, -1>
{
  static bool test()
  {
    return true;
  }
};


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

  auto defaultTests = make_tuple(
    /* 0*/ make_pair( Dune::P0LocalFiniteElement<double,double,2>(simplexGeometry), "p0lfem"),
    /* 1*/ make_pair( Dune::P1LocalFiniteElement<double,double,1>(), "p11dlfem"),
    /* 2*/ make_pair( Dune::P1LocalFiniteElement<double,double,2>(), "p12dlfem"),
    /* 3*/ make_pair( Dune::P1LocalFiniteElement<double,double,3>(), "p13dlfem"),
    /* 4*/ make_pair( Dune::Q1LocalFiniteElement<double,double,1>(), "q11dlfem"),
    /* 5*/ make_pair( Dune::Q1LocalFiniteElement<double,double,2>(), "q12dlfem"),
    /* 6*/ make_pair( Dune::Q1LocalFiniteElement<double,double,3>(), "q13dlfem"),
    /* 7*/ make_pair( Dune::PQ22DLocalFiniteElement<double,double>(simplexGeometry), "pq22dlfem"),
    /* 8*/ make_pair( Dune::RefinedP1LocalFiniteElement<double,double,1>(), "refp11dlfem"),
    /* 9*/ make_pair( Dune::RefinedP1LocalFiniteElement<double,double,2>(), "refp12dlfem"),
    /*10*/ make_pair( Dune::RefinedP1LocalFiniteElement<double,double,3>(), "refp13dlfem"),
    /*11*/ make_pair( Dune::RefinedP0LocalFiniteElement<double,double,1>(), "refp01dlfem"),
    /*12*/ make_pair( Dune::RefinedP0LocalFiniteElement<double,double,2>(), "refp02dlfem"),
    /*13*/ make_pair( Dune::P23DLocalFiniteElement<double,double>(), "p23dlfem"),
    /*14*/ make_pair( Dune::HierarchicalP2LocalFiniteElement<double,double,1>(), "hierarchicalp21dlfem"),
    /*15*/ make_pair( Dune::HierarchicalP2LocalFiniteElement<double,double,2>(), "hierarchicalp22dlfem"),
    /*16*/ make_pair( Dune::HierarchicalP2LocalFiniteElement<double,double,3>(), "hierarchicalp23dlfem"),
    /*17*/ make_pair( Dune::HierarchicalPrismP2LocalFiniteElement<double,double>(), "hierarchicalprismp2lfem"),
    /*18*/ make_pair( Dune::HierarchicalP2WithElementBubbleLocalFiniteElement<double,double,2>(), "hierarchicalp2bubble2dlfem"),
    /*19*/ make_pair( Dune::PrismP1LocalFiniteElement<double,double>(), "prismp1fem"),
    /*20*/ make_pair( Dune::PrismP2LocalFiniteElement<double,double>(), "prismp2fem"),
    /*21*/ make_pair( Dune::QkLocalFiniteElement<double,double,1,1>(), "qk11dlfem"),
    /*22*/ make_pair( Dune::QkLocalFiniteElement<double,double,2,0>(), "qk02dlfem"),
    /*23*/ make_pair( Dune::QkLocalFiniteElement<double,double,2,1>(), "qk12dlfem"),
    /*24*/ make_pair( Dune::QkLocalFiniteElement<double,double,2,2>(), "qk22dlfem"),
    /*25*/ make_pair( Dune::QkLocalFiniteElement<double,double,2,3>(), "qk32dlfem"),
    /*26*/ make_pair( Dune::QkLocalFiniteElement<double,double,3,0>(), "qk03dlfem"),
    /*27*/ make_pair( Dune::QkLocalFiniteElement<double,double,3,1>(), "qk13dlfem"),
    /*28*/ make_pair( Dune::QkLocalFiniteElement<double,double,3,2>(), "qk23dlfem"),
    /*29*/ make_pair( Dune::QkLocalFiniteElement<double,double,3,3>(), "qk33dlfem"),
    /*30*/ make_pair( Dune::DualP1LocalFiniteElement<double,double,1>(), "dualp11dlfem"),
    /*31*/ make_pair( Dune::DualP1LocalFiniteElement<double,double,2>(), "dualp12dlfem"),
    /*32*/ make_pair( Dune::DualP1LocalFiniteElement<double,double,3>(), "dualp13dlfem"),
    /*33*/ make_pair( Dune::DualQ1LocalFiniteElement<double,double,1>(), "dualq11dlfem"),
    /*34*/ make_pair( Dune::DualQ1LocalFiniteElement<double,double,2>(), "dualq12dlfem"),
    /*35*/ make_pair( Dune::DualQ1LocalFiniteElement<double,double,3>(), "dualq13dlfem"),
    /*36*/ make_pair( Dune::RannacherTurekLocalFiniteElement<double,double,2>(), "rannacher_turek2dfem"),
    /*37*/ make_pair( Dune::RannacherTurekLocalFiniteElement<double,double,3>(), "rannacher_turek3dfem"),
    /*38*/ make_pair( Dune::BDM1Cube2DLocalFiniteElement<double,double>(1), "bdm1cube2dlfem"),
    /*39*/ make_pair( Dune::BDM2Cube2DLocalFiniteElement<double,double>(1), "bdm2cube2dlfem"),
    /*40*/ make_pair( Dune::BDM1Simplex2DLocalFiniteElement<double,double>(1), "bdm1simplex2dlfem"),
    /*41*/ make_pair( Dune::RaviartThomasSimplexLocalFiniteElement<2,double,double>(simplexGeometry, 0), "rt0simplex2dlfem"),
    /*42*/ make_pair( Dune::RaviartThomasSimplexLocalFiniteElement<2,double,double>(simplexGeometry, 1), "rt1simplex2dlfem"),
    /*43*/ make_pair( Dune::RaviartThomasCubeLocalFiniteElement<double,double,2,0>(1), "rt0cube2dlfem"),
    /*44*/ make_pair( Dune::RaviartThomasCubeLocalFiniteElement<double,double,3,0>(1), "rt0cube3dlfem"),
    /*45*/ make_pair( Dune::RaviartThomasCubeLocalFiniteElement<double,double,2,1>(1), "rt1cube2dlfem"),
    /*46*/ make_pair( Dune::RaviartThomasCubeLocalFiniteElement<double,double,3,1>(1), "rt1cube3dlfem"),
    /*47*/ make_pair( Dune::RaviartThomasCubeLocalFiniteElement<double,double,2,2>(1), "rt2cube2dlfem"),
    /*48*/ make_pair( Dune::RaviartThomasCubeLocalFiniteElement<double,double,2,3>(1), "rt3cube2dlfem"),
    /*49*/ make_pair( Dune::RaviartThomasCubeLocalFiniteElement<double,double,2,4>(1), "rt4cube2dlfem")
    );

  // Add monomials to the tuple of tests.
  // Monomials are only tested up to order 2 due to the instability of interpolate()
  auto defaultTests_ = Dune::staticAccumulate<0, 2+1>([](auto const k)
  {
      Dune::GeometryType gt;

      return std::make_tuple(
        std::make_pair(Dune::MonomialLocalFiniteElement<double,double,1,k>((gt.makeLine(), gt)),        "monom1d_" + std::to_string(k)),
        std::make_pair(Dune::MonomialLocalFiniteElement<double,double,2,k>((gt.makeTriangle(), gt)),    "monom2d_" + std::to_string(k)),
        std::make_pair(Dune::MonomialLocalFiniteElement<double,double,3,k>((gt.makeTetrahedron(), gt)), "monom3d_" + std::to_string(k))
      );
  }, defaultTests, [](auto&& a, auto&& b) { return std::tuple_cat(std::forward<decltype(a)>(a), std::forward<decltype(b)>(b)); });

  // run tests...
  bool success = Dune::staticAccumulate<0, size(defaultTests_)>([&tests = defaultTests_](auto const i)
  {
    bool b = testFE(std::get<i>(tests).first);
    std::cout << "testFE(" << std::get<i>(tests).second << ") " << (b ? "succeeded\n" : "failed\n");
    return b;
  }, true, [](bool a, bool b) { return a && b; });


  // ---------------------------------------------------------------------------
  // Test some FiniteElements with disabled local-interpolation tests.
  // Intepolation is not yet implemented.

  auto disabledInterpolationTests = make_tuple(
    /* 0*/ make_pair( Dune::BDM1Cube3DLocalFiniteElement<double,double>(1), "bdm1cube3dlfem"),
    /* 1*/ make_pair( Dune::BDM2Simplex2DLocalFiniteElement<double,double>(1), "bdm2simplex2dlfem")
    );

  // run tests...
  success = Dune::staticAccumulate<0, size(disabledInterpolationTests)>([&tests = disabledInterpolationTests](auto const i)
  {
    bool b = testFE(std::get<i>(tests).first, DisableLocalInterpolation);
    std::cout << "testFE(" << std::get<i>(tests).second << ", DisableLocalInterpolation) " << (b ? "succeeded\n" : "failed\n");
    return b;
  }, success, [](bool a, bool b) { return a && b; });


  // ---------------------------------------------------------------------------
  // Test some FiniteElements with disabled jacobian evaluation tests.
  // These pyramid element have only piecewise linear/quadratic basis functions, i.e.
  // not continuousely differentiable in the whole pyramid.

  auto disabledJacobianTests = make_tuple(
    /* 0*/ make_pair( Dune::PyramidP1LocalFiniteElement<double,double>(), "pyramidp1fem"),
    /* 1*/ make_pair( Dune::PyramidP2LocalFiniteElement<double,double>(), "pyramidp2fem")
    );

  // run tests...
  success = Dune::staticAccumulate<0, size(disabledJacobianTests)>([&tests = disabledJacobianTests](auto const i)
  {
    bool b = testFE(std::get<i>(tests).first, DisableJacobian);
    std::cout << "testFE(" << std::get<i>(tests).second << ", DisableJacobian) " << (b ? "succeeded\n" : "failed\n");
    return b;
  }, success, [](bool a, bool b) { return a && b; });

  // ---------------------------------------------------------------------------
  // test various Pk Finite-Elements for variing k

  success = PkLocalFiniteElementTest<1, 1>::test() and success;
  success = PkLocalFiniteElementTest<1, 2>::test() and success;
  success = PkLocalFiniteElementTest<2, 10>::test() and success;
  success = PkLocalFiniteElementTest<3, 10>::test() and success;

  // ---------------------------------------------------------------------------
  // test virtualized FEs

  // run tests...
  success = Dune::staticAccumulate<0, size(defaultTests)>([&tests = defaultTests](auto const i)
  {
    using FE = std::tuple_element_t<0, std::tuple_element_t<i, decltype(defaultTests)>>;
    Dune::LocalFiniteElementVirtualImp<FE> testVirtual(std::get<i>(tests).first);

    bool b1 = testFE(testVirtual);
    std::cout << "testFE(" << std::get<i>(tests).second << "-virtual1) " << (b1 ? "succeeded\n" : "failed\n");

    // notice the second level of virtualization
    Dune::LocalFiniteElementVirtualImp<Dune::LocalFiniteElementVirtualImp<FE>> testVirtual2(testVirtual);

    bool b2 = testFE(testVirtual2);
    std::cout << "testFE(" << std::get<i>(tests).second << "-virtual2) " << (b2 ? "succeeded\n" : "failed\n");

    using Interface = Dune::LocalFiniteElementVirtualInterface<typename FE::Traits::LocalBasisType::Traits>;
    bool b3 = testFE(static_cast<const Interface&>(testVirtual));
    std::cout << "testFE(" << std::get<i>(tests).second << "-virtual3) " << (b3 ? "succeeded\n" : "failed\n");

    return b1 && b2 && b3;
  }, success, [](bool a, bool b) { return a && b; });


  return success ? 0 : 1;
}
catch (Dune::Exception e)
{
  std::cout << e << std::endl;
  return 1;
}
