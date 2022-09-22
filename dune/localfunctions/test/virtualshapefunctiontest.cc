// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#include <config.h>

#undef DUNE_VIRTUAL_SHAPEFUNCTIONS

#include <array>
#include <cstddef>
#include <iostream>

#include <dune/geometry/type.hh>
#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>

#include <dune/localfunctions/lagrange/p0.hh>
#include <dune/localfunctions/lagrange/lagrangesimplex.hh>
#include <dune/localfunctions/lagrange/pq22d.hh>
#include <dune/localfunctions/monomial.hh>

/** \file
    \brief Test the dynamically polymorphic shape function interface

    This file mainly tests whether the polymorphic interface can be properly
    instantiated, compiled and run without crashed.  It does _not_ test whether
    the shape function sets behave correctly.
 */

using namespace Dune;

template <class T>
void syntax_check( const T& )
{}

// A test function to test the local interpolation
template <class D, class R>
struct TestFunction
{
  using DomainType = D;
  using RangeType = R;

  struct Traits {
    using DomainType = D;
    using RangeType = R;
  };

  void evaluate(const DomainType& in, RangeType& out) const {
    // May not be flexible enough to compile for all range types
    out = 1;
  }
};


template <class T>
void testLocalBasis(const LocalBasisVirtualInterface<T>* localBasis)
{
  // call each method once to test that it's there
  syntax_check<unsigned int>( localBasis->order() );
  [[maybe_unused]] unsigned int size = localBasis->size();

  // evaluate the local basis at (0,...,0)
  typename T::DomainType in(0);
  std::vector<typename T::RangeType> out;
  localBasis->evaluateFunction(in, out);
  assert(out.size() == size);

  std::vector<typename T::JacobianType> jacobianOut;
  localBasis->evaluateJacobian(in, jacobianOut);
  assert(jacobianOut.size() == localBasis->size());
}

void testLocalCoefficients(const LocalCoefficientsVirtualInterface* localCoefficients)
{
  if (!localCoefficients)
    DUNE_THROW(Dune::Exception, "Received an invalid pointer to LocalCoefficientsVirtualInterface");

  if (localCoefficients->size() < 1)
    DUNE_THROW(Dune::Exception, "LocalCoefficients does not provide any coefficients!");

  for (std::size_t i=0; i<localCoefficients->size(); i++) {

    // Test the localKey method
    // We just test whether the interface is there.  Correctness is tested elsewhere
    syntax_check<unsigned int>( localCoefficients->localKey(i).subEntity() );
    syntax_check<unsigned int>( localCoefficients->localKey(i).codim() );
    syntax_check<unsigned int>( localCoefficients->localKey(i).index() );

  }
}

template <class DomainType, class RangeType>
void testLocalInterpolation(const LocalInterpolationVirtualInterface<DomainType,RangeType>* localInterpolation)
{
  // Test interpolation of a function object derived from VirtualFunction
  TestFunction<DomainType,RangeType> testFunction;
  std::vector<typename RangeType::field_type> coefficients;

  //////////////////////////////////////////////////////////////////////////////
  //  Part A: Feed the function to the 'interpolate' method in form of
  //    a class providing an evaluate() method.
  //    This way is deprecated since dune-localfunctions 2.7.
  //////////////////////////////////////////////////////////////////////////////
  localInterpolation->interpolate(testFunction, coefficients);

  //////////////////////////////////////////////////////////////////////////////
  //  Part B: Redo the same test, but feed the function to the
  //    'interpolate' method in form of a callable.
  //////////////////////////////////////////////////////////////////////////////
  auto callableTestFunction = [&](const auto& x) {
    RangeType y(0);
    testFunction.evaluate(x,y);
    return y;
  };
  localInterpolation->interpolate(callableTestFunction, coefficients);
}

// Test all methods of a local finite element given as a pointer to the abstract base class
template <class T>
void testLocalFiniteElement(const LocalFiniteElementVirtualInterface<T>* localFiniteElement)
{
  // Test method type()
  std::cout << "Testing local finite element for a " << localFiniteElement->type() << "." << std::endl;

  typedef LocalFiniteElementVirtualInterface<T> FEType;

  // Test the local basis
  const typename FEType::Traits::LocalBasisType& basis = localFiniteElement->localBasis();
  testLocalBasis(&basis);

  // Test the local coefficients
  const typename FEType::Traits::LocalCoefficientsType& coeffs = localFiniteElement->localCoefficients();
  testLocalCoefficients(&coeffs);

  // Test the interpolation
  const typename FEType::Traits::LocalInterpolationType& interp = localFiniteElement->localInterpolation();
  testLocalInterpolation(&interp);

  // Test cloning
  const LocalFiniteElementVirtualInterface<T>* other = localFiniteElement->clone();

  // Make sure new object has the same type
  if (typeid(other).hash_code() != typeid(localFiniteElement).hash_code())
    DUNE_THROW(Dune::Exception, "'clone' method returns object of wrong type");

  // But it shouldn't be the same object
  if (other==localFiniteElement)
    DUNE_THROW(Dune::Exception, "'clone' method returned the object it was called for!");

  delete other;

}

int main (int argc, char *argv[])
{

  typedef Dune::LagrangeSimplexLocalFiniteElement<double, double, 2, 1>::Traits::LocalBasisType::Traits LBTraits;

  const Dune::P0LocalFiniteElement<double, double, 2> p0FE(Dune::GeometryTypes::cube(2));
  const Dune::LocalFiniteElementVirtualImp<Dune::P0LocalFiniteElement<double, double, 2> > p0VFE(p0FE);
  testLocalFiniteElement<LBTraits>(&p0VFE);

  const Dune::PQ22DLocalFiniteElement<double, double> pq2FE(Dune::GeometryTypes::cube(2));
  const Dune::PQ22DLocalFiniteElement<double, double> pq2FE2(pq2FE);

  const Dune::LocalFiniteElementVirtualImp<Dune::PQ22DLocalFiniteElement<double, double> > pq2VFE(pq2FE);
  testLocalFiniteElement<LBTraits>(&pq2VFE);

  const Dune::LocalFiniteElementVirtualImp<Dune::LagrangeSimplexLocalFiniteElement<double, double, 2, 1>> p1VFE;
  testLocalFiniteElement<LBTraits>(&p1VFE);

  const Dune::LocalFiniteElementVirtualImp<Dune::LocalFiniteElementVirtualImp<
    Dune::LagrangeSimplexLocalFiniteElement<double, double, 2, 1>>> p1VVFE;
  testLocalFiniteElement<LBTraits>(&p1VVFE);

  typedef Dune::MonomialLocalFiniteElement<double, double, 2, 7> Monom7;
  const Monom7 monom7FE(Dune::GeometryTypes::cube(2));
  const Dune::LocalFiniteElementVirtualImp<Monom7> monom7VFE(monom7FE);
  const Dune::LocalFiniteElementVirtualImp<
      Dune::LocalFiniteElementVirtualImp<Monom7> > monom7VVFE(monom7VFE);
  testLocalFiniteElement<LBTraits>(&monom7VVFE);

  return 0;

}
