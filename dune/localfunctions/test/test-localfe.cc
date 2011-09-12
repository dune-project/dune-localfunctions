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
#include "../lagrange/q22d.hh"

#include "../refined/refinedp1.hh"
#include "../refined/refinedp0.hh"
#include "../hierarchical/hierarchicalp2.hh"
#include "../hierarchical/hierarchicalp2withelementbubble.hh"
#include "../hierarchical/hierarchicalprismp2.hh"
#include "../rannacherturek/rannacherturek2d.hh"
#include "../raviartthomas/raviartthomas02d.hh"
#include "../raviartthomas/raviartthomas0q2d.hh"
#include "../raviartthomas/raviartthomas0q3d.hh"
#include "../monom.hh"

#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>

#include "test-localfe.hh"

// tmp for testing arbitrary order finite elements
template<unsigned int d, int k>
struct PkLocalFiniteElementTest
{
  static bool test()
  {
    bool success = true;

    Dune::PkLocalFiniteElement<double,double,d,k> pklfem;
    success = testFE(pklfem) and success;

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

template<int k>
bool testMonomials()
{
  bool success = true;
  Dune::GeometryType gt;

  gt.makeLine();
  Dune::MonomLocalFiniteElement<double,double,1,k> monom1d(gt);
  success = testFE(monom1d) and success;

  gt.makeTriangle();
  Dune::MonomLocalFiniteElement<double,double,2,k> monom2d(gt);
  success = testFE(monom2d) and success;

  gt.makeTetrahedron();
  Dune::MonomLocalFiniteElement<double,double,3,k> monom3d(gt);
  success = testFE(monom3d) and success;

  return testMonomials<k-1>() and success;
}

template<>
bool testMonomials<-1>()
{
  return true;
}

int main(int argc, char** argv) try
{
  bool success = true;

  Dune::P0LocalFiniteElement<double,double,2> p0lfem(Dune::GeometryType::simplex);
  success = testFE(p0lfem) and success;

  Dune::P1LocalFiniteElement<double,double,1> p11dlfem;
  success = testFE(p11dlfem) and success;

  Dune::P1LocalFiniteElement<double,double,2> p12dlfem;
  success = testFE(p12dlfem) and success;

  Dune::P1LocalFiniteElement<double,double,3> p13dlfem;
  success = testFE(p13dlfem) and success;

  Dune::Q1LocalFiniteElement<double,double,1> q11dlfem;
  success = testFE(q11dlfem) and success;

  Dune::Q1LocalFiniteElement<double,double,2> q12dlfem;
  success = testFE(q12dlfem) and success;

  Dune::Q1LocalFiniteElement<double,double,3> q13dlfem;
  success = testFE(q13dlfem) and success;

  Dune::PQ22DLocalFiniteElement<double,double> pq22dlfem(
    Dune::GeometryType( Dune::GeometryType::simplex,2) );
  success = testFE(pq22dlfem) and success;

  Dune::RefinedP1LocalFiniteElement<double,double,2> refp12dlfem;
  success = testFE(refp12dlfem) and success;

  Dune::RefinedP1LocalFiniteElement<double,double,3> refp13dlfem;
  success = testFE(refp13dlfem) and success;

  Dune::RefinedP0LocalFiniteElement<double,double,2> refp02dlfem;
  success = testFE(refp02dlfem) and success;

  Dune::P23DLocalFiniteElement<double,double> p23dlfem;
  success = testFE(p23dlfem) and success;

  //    Dune::HierarchicalP2LocalFiniteElement<double,double,1> hierarchicalp21dlfem;
  //    success = testFE(hierarchicalp21dlfem) and success;

  Dune::HierarchicalP2LocalFiniteElement<double,double,2> hierarchicalp22dlfem;
  success = testFE(hierarchicalp22dlfem) and success;

  Dune::HierarchicalP2LocalFiniteElement<double,double,3> hierarchicalp23dlfem;
  success = testFE(hierarchicalp23dlfem) and success;

  Dune::HierarchicalPrismP2LocalFiniteElement<double,double> hierarchicalprismp2lfem;
  success = testFE(hierarchicalprismp2lfem) and success;

  Dune::HierarchicalP2WithElementBubbleLocalFiniteElement<double,double,2> hierarchicalp2bubble2dlfem;
  success = testFE(hierarchicalp2bubble2dlfem) and success;

  Dune::PrismP1LocalFiniteElement<double,double> prismp1fem;
  success = testFE(prismp1fem) and success;

  Dune::PrismP2LocalFiniteElement<double,double> prismp2fem;
  success = testFE(prismp2fem) and success;

  Dune::PyramidP1LocalFiniteElement<double,double> pyramidp1fem;
  success = testFE(pyramidp1fem) and success;

  Dune::PyramidP2LocalFiniteElement<double,double> pyramidp2fem;
  success = testFE(pyramidp2fem) and success;

  success = PkLocalFiniteElementTest<1, 1>::test() and success;

  success = PkLocalFiniteElementTest<2, 10>::test() and success;

  success = PkLocalFiniteElementTest<3, 10>::test() and success;

  Dune::Q22DLocalFiniteElement<double,double> q22dlfem;
  success = testFE(q22dlfem) and success;

  Dune::RT02DLocalFiniteElement<double,double> rt02dlfem(1);
  success = testFE(rt02dlfem) and success;

  Dune::RT0Q2DLocalFiniteElement<double,double> rt0q2dlfem(1);
  success = testFE(rt0q2dlfem) and success;

  Dune::RT0Q3DLocalFiniteElement<double,double> rt0q3dlfem(1);
  success = testFE(rt0q3dlfem) and success;

  Dune::RannacherTurek2DLocalFiniteElement<double,double> rannacher_turek2dfem;
  success = testFE(rannacher_turek2dfem) and success;

  std::cout << "Monomials are only tested up to order 2 due to the instability of interpolate()." << std::endl;
  success = testMonomials<2>() and success;

  // test virtualalized FEs
  // notice that testFE add another level of virtualization
  Dune::LocalFiniteElementVirtualImp< Dune::P1LocalFiniteElement<double,double, 2> >
  p12dlfemVirtual(p12dlfem);
  success = testFE(p12dlfemVirtual) and success;

  Dune::LocalFiniteElementVirtualImp< Dune::PQ22DLocalFiniteElement<double,double> >
  pq22dlfemVirtual(pq22dlfem);
  success = testFE(pq22dlfemVirtual) and success;

  Dune::LocalFiniteElementVirtualImp<
      Dune::LocalFiniteElementVirtualImp<
          Dune::P1LocalFiniteElement<double,double, 2> > >
  p12dlfemVirtualVirtual(p12dlfemVirtual);
  success = testFE(p12dlfemVirtualVirtual) and success;

  Dune::LocalFiniteElementVirtualImp<
      Dune::LocalFiniteElementVirtualImp<
          Dune::PQ22DLocalFiniteElement<double,double> > >
  pq22dlfemVirtualVirtual(pq22dlfemVirtual);
  success = testFE(pq22dlfemVirtualVirtual) and success;

  typedef Dune::LocalFiniteElementVirtualInterface< Dune::P1LocalFiniteElement<double,double, 2>::Traits::LocalBasisType::Traits > Interface;
  success = testFE<Interface>(p12dlfemVirtual) and success;

  return success ? 0 : 1;
}
catch (Dune::Exception e)
{
  std::cout << e << std::endl;
  return 1;
}
