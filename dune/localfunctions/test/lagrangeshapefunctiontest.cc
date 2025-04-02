// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#include <iostream>
#include <typeinfo>
#include <fenv.h>

#include <dune/common/classname.hh>
#include <dune/common/deprecated.hh>
#include <dune/common/test/testsuite.hh>

#include <dune/localfunctions/lagrange/p0.hh>
#include <dune/localfunctions/lagrange/lagrangelfecache.hh>
#include <dune/localfunctions/lagrange/lagrangecube.hh>
#include <dune/localfunctions/lagrange/lagrangeprism.hh>
#include <dune/localfunctions/lagrange/lagrangepyramid.hh>
#include <dune/localfunctions/lagrange/lagrangesimplex.hh>

#define DUNE_DISABLE_DEPRECATION_WARNING_PQ22D
#include <dune/localfunctions/lagrange/pq22d.hh>
#undef DUNE_DISABLE_DEPRECATION_WARNING_PQ22D

#include <dune/localfunctions/test/test-localfe.hh>

/** \file
    \brief Performs some tests for the Lagrange shape functions
 */

double epsilon = 1e-14;
double sqrt_epsilon = std::sqrt(epsilon);

using namespace Dune;

// Generate list of Lagrange points for dim-dimensional simplex
template <int dim>
void getPkTestPoints(unsigned order, unsigned level, std::vector<FieldVector<double,dim> >& test_points)
{
  for (unsigned i = 0; i <= order - level; ++i)
  {
    std::vector<FieldVector<double,dim-1> > test_points_lower_dim;
    getPkTestPoints(order, level + i, test_points_lower_dim);
    double coord = double(i) / order;
    for (unsigned j = 0; j < test_points_lower_dim.size(); ++j)
    {
      FieldVector<double,dim> pos;
      for (int k = 0; k < dim-1; ++k)
        pos[k] = test_points_lower_dim[j][k];
      pos[dim-1] = coord;
      test_points.push_back(pos);
    }
  }
}

// Template specialization to terminate recursion
template <>
void getPkTestPoints(unsigned order, unsigned level, std::vector<FieldVector<double,0> >& test_points)
{
  FieldVector<double,0> pos;
  test_points.push_back(pos);
}

template <class FE>
bool testPk(const FE& local_fe)
{
  const int dim = FE::Traits::LocalBasisType::Traits::dimDomain;
  const unsigned order = local_fe.localBasis().order();

  std::vector<FieldVector<double,1> > values;

  std::vector<FieldVector<double,dim> > test_points;
  getPkTestPoints(order, 0, test_points);

  for (unsigned n = 0; n < test_points.size(); ++n)
  {
    FieldVector<double,dim> pos = test_points[n];

    //////////////////////////////////////////////////////////////////
    //  Verify that shape functions fulfill \phi_i(x_j) = \delta_{ij}
    //  We assume that the shape functions are ordered corresponding
    //  to the test points returned by getPkTestPoints()
    //////////////////////////////////////////////////////////////////

    local_fe.localBasis().evaluateFunction(pos, values);
    for (unsigned i = 0; i < values.size(); ++i)
      if (std::abs(values[i] - double(i==n)) > epsilon)
      {
        std::cerr << "Bug in shape function in local finite element type"
                  << typeid(FE).name() << std::endl;
        std::cerr << "Shape function " << n << " has value " << values[i]
                  << " at position " << pos << " while " << double(i==n)
                  << " was expected" << std::endl;
        return false;
      }

  }

  return true;
}

int main (int argc, char *argv[])
{
#if __linux__ \
  && (!defined __INTEL_COMPILER || __INTEL_COMPILER >= 1010) \
  && (!defined __clang__)
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

  auto testSuite = Dune::TestSuite();

  //////////////////////////////////////////////////////////
  //   Test for the Lagrange property
  //////////////////////////////////////////////////////////

  Dune::Hybrid::forEach(std::index_sequence<1,2,3>{},[&](auto dim)
  {
    Dune::Hybrid::forEach(std::index_sequence<1,2,3,4>{},[&](auto order)
    {
      auto lfe = LagrangeSimplexLocalFiniteElement<double,double,dim,order>();
      testSuite.check(testPk(lfe))
        << "Lagrange property not satisfied for " << Dune::className(lfe);
    });
  });

  //////////////////////////////////////////////////////////
  //   Run the standard tests
  //////////////////////////////////////////////////////////

  // This check extends the generic testFE() by additionally
  // wrapping the LFE into the virtual interface twice
  // Check local finite element and its virtualized variants
  // Notice that testFE add another level of virtualization
  auto testVirtualLFE = [&](const auto& lfe, auto... args)
  {
    auto testSuite = Dune::TestSuite(Dune::className(lfe));
    using LFE = std::decay_t<decltype(lfe)>;
    using LocalBasisTraits = typename LFE::Traits::LocalBasisType::Traits;
    using Interface = Dune::LocalFiniteElementVirtualInterface<LocalBasisTraits>;
    auto vlfe = Dune::LocalFiniteElementVirtualImp<LFE>(lfe);
    auto vvlfe = Dune::LocalFiniteElementVirtualImp<decltype(vlfe)>(vlfe);
    const auto& ilfe = static_cast<const Interface&>(vlfe);
    testSuite.check(testFE(lfe, args...)) << "Check of raw local finite element";
    testSuite.check(testFE(vlfe, args...)) << "Check of virtualized local finite element";
    testSuite.check(testFE(vvlfe, args...)) << "Check of double virtualized local finite element";
    testSuite.check(testFE(ilfe, args...)) << "Check of virtualized local finite element via interface";
    return testSuite;
  };

  // Special implementations

  auto p0LFE = Dune::P0LocalFiniteElement<double,double,2>(GeometryTypes::simplex(2));
  testSuite.subTest(testVirtualLFE(p0LFE, DisableNone, 0));

DUNE_NO_DEPRECATED_BEGIN
  auto mixedPQ22DLFE = Dune::PQ22DLocalFiniteElement<double,double>(GeometryTypes::simplex(2));
  testSuite.subTest(testVirtualLFE(mixedPQ22DLFE, DisableNone, 0));
DUNE_NO_DEPRECATED_END

  auto prismP1LFE = LagrangePrismLocalFiniteElement<double,double,1>();
  testSuite.subTest(testVirtualLFE(prismP1LFE, DisableNone, 2));

  auto prismP2LFE = LagrangePrismLocalFiniteElement<double,double,2>();
  testSuite.subTest(testVirtualLFE(prismP2LFE, DisableNone, 1));

  // Pyramid shapefunctions are not differentiable on the plane where xi[0]=xi[1].
  // So let's skip test points on this plane
  auto xySkip = [](const FieldVector<double,3>& xi){return std::abs(xi[0]-xi[1]) < 1e-8;};

  auto pyramidP1LFE = LagrangePyramidLocalFiniteElement<double,double,1>();
  testSuite.subTest(testVirtualLFE(pyramidP1LFE, DisableNone, 1, xySkip));

  auto pyramidP2LFE = LagrangePyramidLocalFiniteElement<double,double,2>();
  testSuite.subTest(testVirtualLFE(pyramidP2LFE, DisableNone, 1, xySkip));

  // Simplex implementations
  Dune::Hybrid::forEach(std::index_sequence<1,2,3>{},[&](auto dim)
  {
    Dune::Hybrid::forEach(std::make_index_sequence<5>{},[&](auto order)
    {
      auto diffOrder = 2;
      auto lfe = LagrangeSimplexLocalFiniteElement<double,double,dim,order>();
      testSuite.subTest(testVirtualLFE(lfe, DisableNone, diffOrder));
    });
  });

  // Cube implementations
  Dune::Hybrid::forEach(std::index_sequence<1,2,3>{},[&](auto dim)
  {
    Dune::Hybrid::forEach(std::make_index_sequence<5>{},[&](auto order)
    {
      auto diffOrder = 2;
      auto lfe = LagrangeCubeLocalFiniteElement<double,double,dim,order>();
      testSuite.subTest(testVirtualLFE(lfe, DisableNone, diffOrder));
    });
  });

  // Test the LagrangeLocalFiniteElementCache
  auto lagrangeLFECache = LagrangeLocalFiniteElementCache<double,double,2,2>();
  testSuite.subTest(testVirtualLFE(lagrangeLFECache.get(GeometryTypes::simplex(2))));
  testSuite.subTest(testVirtualLFE(lagrangeLFECache.get(GeometryTypes::cube(2))));

  // Test whether asking the cache for an element of the wrong dimension throws an exception
  testSuite.checkThrow([&]{
    lagrangeLFECache.get(GeometryTypes::simplex(1));
  }) << "LagrangeLocalFiniteElementCache does not throw for non-existing element";

  return testSuite.exit();
}
