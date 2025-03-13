// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#include <cstddef>
#include <iostream>
#include <typeinfo>
#include <cstdlib>
#include <vector>

#include <dune/common/gmpfield.hh>

// Lagrange type elements
#include <dune/localfunctions/lagrange.hh>
#include <dune/localfunctions/lagrange/equidistantpoints.hh>

// DG type elements
#include <dune/localfunctions/orthonormal.hh>

// Raviart Thomas type elements
#include <dune/localfunctions/raviartthomas.hh>

#include "test-localfe.hh"



// Check if the traces of the basis functions on all
// faces have the desired polynomial order by comparing
// them to the interpolation into an appropriate trace
// FE space at a set of sample points.
template<class FE>
bool checkTraceOrder(const FE& fe, int order, double tol)
{
  bool success = true;

  using namespace Dune::Indices;
  using T = typename FE::Traits::LocalBasisType::Traits::RangeFieldType;
  using Range = typename FE::Traits::LocalBasisType::Traits::RangeType;
  constexpr auto dim = FE::Traits::LocalBasisType::Traits::dimDomain;

  auto evaluationBuffer = std::vector<Range>();
  auto coefficients = std::vector<Range>();

  // Loop over all faces with 1<=codim<dim
  const auto& re = Dune::referenceElement<T, dim>(fe.type());
  Dune::Hybrid::forEach(Dune::range(_1, Dune::index_constant<dim>{}), [&](auto codim) {
    for(auto face : Dune::range(re.size(codim)))
    {
      // Check if the trace on the face is in the appropriate
      // polynomial space by interpolating each basis function
      // into a corresponding polynomial trace space and comparing
      // the results.
      using TraceFE = typename Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet, dim-codim, T, T>;
      using std::fabs;
      auto faceGeometry = re.template geometry<codim>(face);
      auto traceFE = TraceFE(re.type(face, codim), order);
      const auto& sampleRule = Dune::QuadratureRules<T, dim-codim>::rule(re.type(face, codim), 2*order);

      for(auto i : Dune::range(fe.size()))
      {
        auto fe_i = [&](auto x) {
          fe.localBasis().evaluateFunction(faceGeometry.global(x), evaluationBuffer);
          return evaluationBuffer[i];
        };

        traceFE.localInterpolation().interpolate(fe_i, coefficients);
        auto fe_i_interpolated = [&](auto x) {
          traceFE.localBasis().evaluateFunction(x, evaluationBuffer);
          auto y = Range(0);
          for(auto j : Dune::range(traceFE.size()))
            y += evaluationBuffer[j]*coefficients[j];
          return y;
        };
        double mismatch = 0;
        for(auto [x, dummy] : sampleRule)
          mismatch = std::max(mismatch, double(fabs(fe_i(x) - fe_i_interpolated(x))));
        if (mismatch > tol)
        {
          std::cout
            << "The trace of basis function " << i
            << " on face " << face
            << " with codim " << codim
            << " of " << fe.type()
            << " is not in the desired polynomial space."
            << " Maximal mismatch to interpolation is " << mismatch
            << std::endl;
          success = false;
        }
      }
    }
  });
  return success;
}



int main(int argc, char** argv) try
{
  bool success = true;

  {
    double tol = 1e-7;
    using FE2D = Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet, 2, double, double>;
    using FE3D = Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet, 3, double, double>;
    for(auto order : {1,2,3,4})
    {
      std::cout << "Checking traces for Lagrange bases of order " << order << std::endl;
      success &= checkTraceOrder(FE2D(Dune::GeometryTypes::triangle, order), order, tol);
      success &= checkTraceOrder(FE2D(Dune::GeometryTypes::quadrilateral, order), order, tol);
      success &= checkTraceOrder(FE3D(Dune::GeometryTypes::tetrahedron, order), order, tol);
      success &= checkTraceOrder(FE3D(Dune::GeometryTypes::hexahedron, order), order, tol);
      success &= checkTraceOrder(FE3D(Dune::GeometryTypes::prism, order), order, tol);
      success &= checkTraceOrder(FE3D(Dune::GeometryTypes::pyramid, order), order, tol);
    }
  }

  std::cout << "Testing LagrangeLocalFiniteElement<EquidistantPointSet> on 3d"
            << " simplex elements with double precision" << std::endl;
  for (unsigned int order : {1, 2, 4, 6})
  {
    std::cout << "order : " << order << std::endl;
    Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet,3,double,double>
    lagrangeSimplex(Dune::GeometryTypes::simplex(3), order);
    TEST_FE(lagrangeSimplex);
  }
  std::cout << "Testing LagrangeLocalFiniteElement<EquidistantPointSet> on 2d"
            << " cube elements with double precision" << std::endl;
  for (unsigned int order : {1, 2, 4})
  {
    std::cout << "order : " << order << std::endl;
    Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet,2,double,double>
    lagrangeCube(Dune::GeometryTypes::cube(2), order);
    TEST_FE(lagrangeCube);
  }
#if HAVE_GMP
  std::cout << "Testing LagrangeLocalFiniteElement<EquidistantPointSet> on 2d"
            << " simplex elements with higher precision" << std::endl;
  for (unsigned int order : {5, 8})
  {
    std::cout << "order : " << order << std::endl;
    Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet,2,double,double,
        Dune::GMPField<64>,Dune::GMPField<256> >
    lagrangeSimplex(Dune::GeometryTypes::simplex(2), order);
    TEST_FE(lagrangeSimplex);
  }
#endif
  std::cout << "Testing DGLagrangeLocalFiniteElement<EquidistantPointSet> on 3d"
            << " cube elements with double precision" << std::endl;
  for (unsigned int order : {1, 2})
  {
    std::cout << "order : " << order << std::endl;
    typedef Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet,3,double,double> FE;
    Dune::DGLocalFiniteElement<FE> dglagrangeCube(Dune::GeometryTypes::cube(3), order);
    TEST_FE(dglagrangeCube);
  }
  std::cout << "Testing L2LagrangeLocalFiniteElement<EquidistantPointSet> on 3d"
            << " cube elements with double precision" << std::endl;
  for (unsigned int order : {2, 3})
  {
    std::cout << "order : " << order << std::endl;
    typedef Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet,3,double,double> FE;
    Dune::L2LocalFiniteElement<FE> dglagrangeCube(Dune::GeometryTypes::cube(3), order);
    TEST_FE(dglagrangeCube);
  }
#if HAVE_GMP
  std::cout << "Testing OrthonormalFiniteElement on 3d"
            << " prism elements with higher precision" << std::endl;
  for (unsigned int order : {6})
  {
    std::cout << "order : " << order << std::endl;
    Dune::OrthonormalLocalFiniteElement<3,double,double,
        Dune::GMPField<64>,Dune::GMPField<256> >
    onbPrism(Dune::GeometryTypes::prism, order);
    TEST_FE(onbPrism);
  }
#endif
  std::cout << "Testing OrthonormalFiniteElement on 3d"
            << " prism elements with double precision" << std::endl;
  for (unsigned int order : {4, 2})
  {
    std::cout << "order : " << order << std::endl;
    Dune::OrthonormalLocalFiniteElement<3,double,double>
    onbPrism(Dune::GeometryTypes::prism, order);
    TEST_FE(onbPrism);
  }
  std::cout << "Testing RaviartThomasSimplexFiniteElement on 3d"
            << " simplex elements with double precision" << std::endl;
  for (unsigned int order : {0, 1, 4})
  {
    std::cout << "order : " << order << std::endl;
    Dune::RaviartThomasSimplexLocalFiniteElement<3,double,double>
    rtSimplex(Dune::GeometryTypes::simplex(3), order);
    TEST_FE(rtSimplex);
  }
  return success ? 0 : 1;
}
catch (const Dune::Exception &e)
{
  std::cout << e << std::endl;
  return 1;
}
