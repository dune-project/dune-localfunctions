// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstddef>
#include <iostream>
#include <typeinfo>
#include <cstdlib>
#include <vector>

#include <dune/common/function.hh>

#include <dune/common/gmpfield.hh>

// Lagrange type elements
#include <dune/localfunctions/lagrange.hh>
#include <dune/localfunctions/lagrange/equidistantpoints.hh>

// DG type elements
#include <dune/localfunctions/orthonormal.hh>

// Raviart Thomas type elements
#include <dune/localfunctions/raviartthomas.hh>

#include "test-localfe.hh"

int main(int argc, char** argv) try
{
  bool success = true;

  std::cout << "Testing LagrangeLocalFiniteElement<EquidistantPointSet> on 3d"
            << " simplex elements with double precision" << std::endl;
  for (unsigned int order=1; order<=6; ++order)
  {
    std::cout << "order : " << order << std::endl;
    Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet,3,double,double>
    lagrangeSimplex(Dune::GeometryType(Dune::GeometryType::simplex, 3), order);
    TEST_FE(lagrangeSimplex);
  }
  std::cout << "Testing LagrangeLocalFiniteElement<EquidistantPointSet> on 2d"
            << " cube elements with double precision" << std::endl;
  for (unsigned int order=1; order<=4; ++order)
  {
    std::cout << "order : " << order << std::endl;
    Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet,2,double,double>
    lagrangeCube(Dune::GeometryType(Dune::GeometryType::cube, 2), order);
    TEST_FE(lagrangeCube);
  }
#if HAVE_GMP
  std::cout << "Testing LagrangeLocalFiniteElement<EquidistantPointSet> on 2d"
            << " simplex elements with higher precision" << std::endl;
  for (unsigned int order=5; order<=8; ++order)
  {
    std::cout << "order : " << order << std::endl;
    Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet,2,double,double,
        Dune::GMPField<64>,Dune::GMPField<256> >
    lagrangeSimplex(Dune::GeometryType(Dune::GeometryType::simplex, 2), order);
    TEST_FE(lagrangeSimplex);
  }
#endif
  std::cout << "Testing DGLagrangeLocalFiniteElement<EquidistantPointSet> on 3d"
            << " cube elements with double precision" << std::endl;
  for (unsigned int order=1; order<=2; ++order)
  {
    std::cout << "order : " << order << std::endl;
    typedef Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet,3,double,double> FE;
    Dune::DGLocalFiniteElement<FE> dglagrangeCube(Dune::GeometryType(Dune::GeometryType::cube, 3), order);
    TEST_FE(dglagrangeCube);
  }
  std::cout << "Testing L2LagrangeLocalFiniteElement<EquidistantPointSet> on 3d"
            << " cube elements with double precision" << std::endl;
  for (unsigned int order=1; order<=3; ++order)
  {
    std::cout << "order : " << order << std::endl;
    typedef Dune::LagrangeLocalFiniteElement<Dune::EquidistantPointSet,3,double,double> FE;
    Dune::L2LocalFiniteElement<FE> dglagrangeCube(Dune::GeometryType(Dune::GeometryType::cube, 3), order);
    TEST_FE(dglagrangeCube);
  }
#if HAVE_GMP
  std::cout << "Testing OrthonormalFiniteElement on 3d"
            << " prism elements with higher precision" << std::endl;
  for (unsigned int order=7; order>=4; --order)
  {
    std::cout << "order : " << order << std::endl;
    Dune::OrthonormalLocalFiniteElement<3,double,double,
        Dune::GMPField<64>,Dune::GMPField<256> >
    onbPrism(Dune::GeometryType(Dune::GeometryType::prism, 3), order);
    TEST_FE(onbPrism);
  }
#endif
  std::cout << "Testing OrthonormalFiniteElement on 3d"
            << " prism elements with double precision" << std::endl;
  for (unsigned int order=4; order>=1; --order)
  {
    std::cout << "order : " << order << std::endl;
    Dune::OrthonormalLocalFiniteElement<3,double,double>
    onbPrism(Dune::GeometryType(Dune::GeometryType::prism, 3), order);
    TEST_FE(onbPrism);
  }
  std::cout << "Testing RaviartThomasSimplexFiniteElement on 3d"
            << " simplex elements with double precision" << std::endl;
  for (unsigned int order=0; order<=4; ++order)
  {
    std::cout << "order : " << order << std::endl;
    Dune::RaviartThomasSimplexLocalFiniteElement<3,double,double>
    rtSimplex(Dune::GeometryType(Dune::GeometryType::simplex, 3), order);
    TEST_FE(rtSimplex);
  }
  return success ? 0 : 1;
}
catch (const Dune::Exception &e)
{
  std::cout << e << std::endl;
  return 1;
}
