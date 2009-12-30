// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <cstddef>
#include <iostream>

#include <dune/localfunctions/common/virtualinterface.hh>

#include <dune/localfunctions/p1.hh>

/** \file
    \brief Test the dynamically polymorphic shape function interface

    This file mainly tests whether the polymorphic interface can be properly
    instantiated, compiled and run without crashed.  It does _not_ test whether
    the shape function sets behave correctly.
 */

using namespace Dune;

// A test function to test the local interpolation
template <class DomainType, class RangeType>
struct TestFunction
  : public VirtualFunction<DomainType,RangeType>
{
  void evaluate(const DomainType& in, RangeType& out) const {
    // May not be flexible enough to compile for all range types
    out = 1;
  }
};


template <class LocalBasisTraits>
void testLocalBasis(const C0LocalBasisVirtualInterface<LocalBasisTraits>* localBasis)
{
  // call each method once to test that it's there
  unsigned int size = localBasis->size();
  unsigned int order = localBasis->order();

  // evaluate the local basis at (0,...,0)
  typename LocalBasisTraits::DomainType in(0);
  std::vector<typename LocalBasisTraits::RangeType> out;
  localBasis->evaluateFunction(in, out);
  assert(out.size() == size);

  // test whether this local basis has a Jacobian and evaluate it if possible
  /** \todo This dynamic testing is Augenwischerei.  If localBasis was not actually
      derived from C1LocalBasisVirtualInterface then LocalBasisTraits would not
      contain a JacobianType and the whole thing wouldn't compile...
   */
  if (dynamic_cast<const C1LocalBasisVirtualInterface<LocalBasisTraits>*>(localBasis)) {

    std::vector<typename LocalBasisTraits::JacobianType> jacobianOut;
    dynamic_cast<const C1LocalBasisVirtualInterface<LocalBasisTraits>*>(localBasis)->evaluateJacobian(in, jacobianOut);
    assert(jacobianOut.size() == size);

  }

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
    unsigned int subEntity = localCoefficients->localKey(i).subEntity();
    unsigned int codim     = localCoefficients->localKey(i).codim();
    unsigned int index     = localCoefficients->localKey(i).index();

  }
}

template <class DomainType, class RangeType>
void testLocalInterpolation(const LocalInterpolationVirtualInterface<DomainType,RangeType>* localInterpolation)
{
  // Test interpolation of a function object derived from VirtualFunction
  TestFunction<DomainType,RangeType> testFunction;
  std::vector<typename RangeType::field_type> coefficients;
  localInterpolation->interpolate(testFunction, coefficients);
}

// Test all methods of a local finite element given as a pointer to the abstract base class
template <class LocalBasisTraits>
void testLocalFiniteElement(const LocalFiniteElementVirtualInterface<LocalBasisTraits>* localFiniteElement)
{
  // Test method type()
  std::cout << "Testing local finite element for a " << localFiniteElement->type() << "." << std::endl;

  // Test the local basis
  const typename LocalFiniteElementVirtualInterface<LocalBasisTraits>::Traits::LocalBasisType& basis = localFiniteElement->localBasis();
  testLocalBasis(&basis);

  // Test the local coefficients
  const typename LocalFiniteElementVirtualInterface<LocalBasisTraits>::Traits::LocalCoefficientsType& coeffs = localFiniteElement->localCoefficients();
  testLocalCoefficients(&coeffs);

  // Test the interpolation
  const typename LocalFiniteElementVirtualInterface<LocalBasisTraits>::Traits::LocalInterpolationType& interp = localFiniteElement->localInterpolation();
  testLocalInterpolation(&interp);

}


int main (int argc, char *argv[]) try
{

  const Dune::LocalFiniteElementVirtualImp<Dune::P1LocalFiniteElement<double, double, 2> > virtualLocalSimplexFE;

  testLocalFiniteElement(&virtualLocalSimplexFE);

  return 0;

}
catch (Exception e) {

  std::cout << e << std::endl;
  return 1;
}
