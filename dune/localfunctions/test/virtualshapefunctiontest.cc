// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <cstddef>
#include <iostream>

#include <dune/localfunctions/common/virtualinterface.hh>

#include <dune/localfunctions/p1.hh>

/** \file
    \brief Test the dynamically polymorphic shape function interface
 */

using namespace Dune;

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

// Test all methods of a local finite element given as a pointer to the abstract base class
template <class LocalBasisTraits>
void testLocalFiniteElement(const LocalFiniteElementVirtualInterface<LocalBasisTraits>* localFiniteElement)
{
  // Test method type()
  std::cout << "Testing local finite element for a " << localFiniteElement->type() << "." << std::endl;

  // Test the local basis

  // Test the local coefficients
  testLocalCoefficients(&localFiniteElement->localCoefficients());

  // Test the interpolation
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
