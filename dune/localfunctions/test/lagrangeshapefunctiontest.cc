// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#include <config.h>

#include <iostream>
#include <typeinfo>
#include <fenv.h>

#include <dune/localfunctions/lagrange/p0.hh>
#include <dune/localfunctions/lagrange/pq22d.hh>
#include <dune/localfunctions/lagrange/lagrangelfecache.hh>
#include <dune/localfunctions/lagrange/lagrangecube.hh>
#include <dune/localfunctions/lagrange/lagrangeprism.hh>
#include <dune/localfunctions/lagrange/lagrangepyramid.hh>
#include <dune/localfunctions/lagrange/lagrangesimplex.hh>

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
    //  Verfiy that shape functions fulfill \phi_i(x_j) = \delta_{ij}
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

  bool success = true;

  //////////////////////////////////////////////////////////
  //   Test for the Lagrange property
  //////////////////////////////////////////////////////////
  LagrangeSimplexLocalFiniteElement<double,double,1,1> p11d;
  success &= testPk(p11d);

  LagrangeSimplexLocalFiniteElement<double,double,2,1> p12d;
  success &= testPk(p12d);

  LagrangeSimplexLocalFiniteElement<double,double,3,1> p13d;
  success &= testPk(p13d);

  //     P23DLocalFiniteElement does not fulfill above assumption on the
  //     ordering of the shape functions
  //     P23DLocalFiniteElement<double,double> p23d;
  //     testPk(p23d);

  LagrangeSimplexLocalFiniteElement<double,double,2,1> pk12d;
  success &= testPk(pk12d);
  LagrangeSimplexLocalFiniteElement<double,double,2,2> pk22d;
  success &= testPk(pk22d);
  LagrangeSimplexLocalFiniteElement<double,double,2,3> pk32d;
  success &= testPk(pk32d);
  LagrangeSimplexLocalFiniteElement<double,double,2,4> pk42d;
  success &= testPk(pk42d);

  LagrangeSimplexLocalFiniteElement<double,double,3,1> pk13d;
  success &= testPk(pk13d);
  LagrangeSimplexLocalFiniteElement<double,double,3,2> pk23d;
  success &= testPk(pk23d);
  LagrangeSimplexLocalFiniteElement<double,double,3,3> pk33d;
  success &= testPk(pk33d);
  LagrangeSimplexLocalFiniteElement<double,double,3,4> pk43d;
  success &= testPk(pk43d);

  //////////////////////////////////////////////////////////
  //   Run the standard tests
  //////////////////////////////////////////////////////////
  P0LocalFiniteElement<double,double,2> p0lfem(
  GeometryTypes::simplex(2));
  TEST_FE(p0lfem);

  LagrangeSimplexLocalFiniteElement<double,double,1,1> p11dlfem;
  TEST_FE3(p11dlfem,DisableNone,2);

  LagrangeSimplexLocalFiniteElement<double,double,2,1> p12dlfem;
  TEST_FE3(p12dlfem,DisableNone,2);

  LagrangeSimplexLocalFiniteElement<double,double,3,1> p13dlfem;
  TEST_FE3(p13dlfem,DisableNone,2);

  LagrangeSimplexLocalFiniteElement<double,double,1,1> q11dlfem;
  TEST_FE(q11dlfem);

  LagrangeSimplexLocalFiniteElement<double,double,2,1> q12dlfem;
  TEST_FE(q12dlfem);

  LagrangeSimplexLocalFiniteElement<double,double,3,1> q13dlfem;
  TEST_FE(q13dlfem);

  LagrangeSimplexLocalFiniteElement<double,double,2,1> p11dfem;
  TEST_FE(p11dfem);

  PQ22DLocalFiniteElement<double,double> pq22dlfem(
    GeometryTypes::simplex(2));
  TEST_FE(pq22dlfem);

  LagrangeSimplexLocalFiniteElement<double,double,3,2> p23dlfem;
  TEST_FE(p23dlfem);

  LagrangePrismLocalFiniteElement<double,double,1> prismp1fem;
  TEST_FE3(prismp1fem, DisableNone, 2);

  LagrangePrismLocalFiniteElement<double,double,2> prismp2fem;
  TEST_FE3(prismp2fem, DisableNone, 1);

  // Pyramid shapefunctions are not differentiable on the plane where xi[0]=xi[1].
  // So let's skip test points on this plane
  auto xySkip = [](const FieldVector<double,3>& xi){return std::abs(xi[0]-xi[1]) < 1e-8;};

  LagrangePyramidLocalFiniteElement<double,double,1> pyramidp1fem;
  TEST_FE4(pyramidp1fem, DisableNone, 1, xySkip);

  LagrangePyramidLocalFiniteElement<double,double,2> pyramidp2fem;
  TEST_FE4(pyramidp2fem, DisableNone, 1, xySkip);

  Hybrid::forEach(std::make_index_sequence<4>{},[&success](auto i)
  {
    LagrangeSimplexLocalFiniteElement<double,double,1,i> pklfem;
    TEST_FE(pklfem);
  });

  Hybrid::forEach(std::make_index_sequence<5>{},[&success](auto i)
  {
    LagrangeSimplexLocalFiniteElement<double,double,2,i> pklfem;
    TEST_FE3(pklfem,DisableNone,2);
  });

  Hybrid::forEach(std::make_index_sequence<6>{},[&success](auto i)
  {
    LagrangeSimplexLocalFiniteElement<double,double,3,i> pklfem;
    TEST_FE(pklfem);
  });

  // --------------------------------------------------------
  //  Test some instantiations of LagrangeCubeLocalFiniteElement
  // --------------------------------------------------------
  LagrangeCubeLocalFiniteElement<double,double,1,1> qk11dlfem;
  TEST_FE3(qk11dlfem,DisableNone,2);

  LagrangeCubeLocalFiniteElement<double,double,2,0> qk02dlfem;
  TEST_FE3(qk02dlfem,DisableNone,2);

  LagrangeCubeLocalFiniteElement<double,double,2,1> qk12dlfem;
  TEST_FE3(qk12dlfem,DisableNone,2);

  LagrangeCubeLocalFiniteElement<double,double,2,2> qk22dlfem;
  TEST_FE3(qk22dlfem,DisableNone,2);

  LagrangeCubeLocalFiniteElement<double,double,2,3> qk32dlfem;
  TEST_FE3(qk32dlfem,DisableNone,2);

  LagrangeCubeLocalFiniteElement<double,double,3,0> qk03dlfem;
  TEST_FE3(qk03dlfem,DisableNone,2);

  LagrangeCubeLocalFiniteElement<double,double,3,1> qk13dlfem;
  TEST_FE3(qk13dlfem,DisableNone,2);

  LagrangeCubeLocalFiniteElement<double,double,3,2> qk23dlfem;
  TEST_FE3(qk23dlfem,DisableNone,2);

  LagrangeCubeLocalFiniteElement<double,double,3,3> qk33dlfem;
  TEST_FE3(qk33dlfem,DisableNone,2);

  // test virtualized FEs
  // notice that testFE add another level of virtualization
  LocalFiniteElementVirtualImp< LagrangeSimplexLocalFiniteElement<double,double,2,1> >
  p12dlfemVirtual(p12dlfem);
  TEST_FE(p12dlfemVirtual);

  LocalFiniteElementVirtualImp< PQ22DLocalFiniteElement<double,double> >
  pq22dlfemVirtual(pq22dlfem);
  TEST_FE(pq22dlfemVirtual);

  LocalFiniteElementVirtualImp<
      LocalFiniteElementVirtualImp<
          LagrangeSimplexLocalFiniteElement<double,double,2,1> > >
  p12dlfemVirtualVirtual(p12dlfemVirtual);
  TEST_FE(p12dlfemVirtualVirtual);

  LocalFiniteElementVirtualImp<
      LocalFiniteElementVirtualImp<
          PQ22DLocalFiniteElement<double,double> > >
  pq22dlfemVirtualVirtual(pq22dlfemVirtual);
  TEST_FE(pq22dlfemVirtualVirtual);

  typedef LocalFiniteElementVirtualInterface< LagrangeSimplexLocalFiniteElement<double,double,2,1>::Traits::LocalBasisType::Traits > Interface;
  TEST_FE(static_cast<const Interface&>(p12dlfemVirtual));

  // Test the LagrangeLocalFiniteElementCache
  LagrangeLocalFiniteElementCache<double,double,2,2> lagrangeLFECache;
  TEST_FE(lagrangeLFECache.get(GeometryTypes::simplex(2)));
  TEST_FE(lagrangeLFECache.get(GeometryTypes::cube(2)));

  // Test whether asking the cache for an element of the wrong dimension throws an exception
  bool lagrangeLFESuccess = false;
  try {
    auto doesntExist = lagrangeLFECache.get(GeometryTypes::simplex(1));
  } catch (Dune::Exception& e)
  {
    lagrangeLFESuccess = true;
  }
  success &= lagrangeLFESuccess;

  return success ? 0 : 1;
}
