// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#undef DUNE_VIRTUAL_SHAPEFUNCTIONS

#include <array>
#include <cstddef>
#include <iostream>

#include <dune/geometry/type.hh>
#include <dune/localfunctions/common/virtualinterface.hh>

#include <dune/localfunctions/lagrange/p0.hh>
#include <dune/localfunctions/lagrange/p1.hh>
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
template <class DomainType, class RangeType>
struct TestFunction
//    : public VirtualFunction<DomainType,RangeType>
  : public Function<const DomainType&,RangeType&>
{
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
  unsigned int size DUNE_UNUSED = localBasis->size();

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
  localInterpolation->interpolate(testFunction, coefficients);
}

template <class Interface, int order>
struct EvaluateTest
{
  static void test(const Interface& fe)
  {
    typedef typename Interface::Traits::LocalBasisType::Traits LBTraits;

    std::array<int,order> d;
    for(unsigned int i=0; i<d.size(); ++i)
      d[i] = 0;

    typename LBTraits::DomainType x;
    x = 0;

    typename std::vector<typename LBTraits::RangeType> y1;
    typename std::vector<typename LBTraits::RangeType> y2;

    fe.localBasis().evaluate(d,x,y1);
    fe.localBasis().template evaluate<order>(d,x,y2);

    for(unsigned int i=0; i<d.size(); ++i)
      if (y1[i] != y2[i])
        DUNE_THROW(Dune::Exception, "result of template evaluate<order>() and virtual evaluate() do not coincide");

    EvaluateTest<Interface, order-1>::test(fe);
  }
};

template <class Interface>
struct EvaluateTest<Interface, -1>
{
  static void test(const Interface& fe)
  {}
};


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

  EvaluateTest<FEType, T::diffOrder>::test(*localFiniteElement);
}

int main (int argc, char *argv[]) try
{

  typedef Dune::P1LocalFiniteElement<double, double, 2>::Traits::LocalBasisType::Traits LBTraits;
  typedef Dune::FixedOrderLocalBasisTraits<LBTraits,0>::Traits C0LBTraits;
  typedef Dune::FixedOrderLocalBasisTraits<LBTraits,1>::Traits C1LBTraits;
  typedef Dune::FixedOrderLocalBasisTraits<LBTraits,2>::Traits C2LBTraits;

  const Dune::P0LocalFiniteElement<double, double, 2> p0FE(Dune::GeometryType(Dune::GeometryType::cube, 2));
  const Dune::LocalFiniteElementVirtualImp<Dune::P0LocalFiniteElement<double, double, 2> > p0VFE(p0FE);
  testLocalFiniteElement<C0LBTraits>(&p0VFE);

  const Dune::PQ22DLocalFiniteElement<double, double> pq2FE(Dune::GeometryType(Dune::GeometryType::cube, 2));
  const Dune::PQ22DLocalFiniteElement<double, double> pq2FE2(pq2FE);

  const Dune::LocalFiniteElementVirtualImp<Dune::PQ22DLocalFiniteElement<double, double> > pq2VFE(pq2FE);
  testLocalFiniteElement<C0LBTraits>(&pq2VFE);

  const Dune::LocalFiniteElementVirtualImp<Dune::P1LocalFiniteElement<double, double, 2> > p1VFE;
  testLocalFiniteElement<C0LBTraits>(&p1VFE);
  testLocalFiniteElement<C1LBTraits>(&p1VFE);
  testLocalFiniteElement<C2LBTraits>(&p1VFE);

  const Dune::LocalFiniteElementVirtualImp<
      Dune::LocalFiniteElementVirtualImp<Dune::P1LocalFiniteElement<double, double, 2> > > p1VVFE;
  testLocalFiniteElement<C0LBTraits>(&p1VVFE);
  testLocalFiniteElement<C1LBTraits>(&p1VVFE);
  testLocalFiniteElement<C2LBTraits>(&p1VVFE);

  typedef Dune::MonomialLocalFiniteElement<double, double, 2, 7> Monom7;
  const Monom7 monom7FE(Dune::GeometryType(Dune::GeometryType::cube, 2));
  const Dune::LocalFiniteElementVirtualImp<Monom7> monom7VFE(monom7FE);
  const Dune::LocalFiniteElementVirtualImp<
      Dune::LocalFiniteElementVirtualImp<Monom7> > monom7VVFE(monom7VFE);
  testLocalFiniteElement<C0LBTraits>(&monom7VVFE);
  testLocalFiniteElement<C1LBTraits>(&monom7VVFE);
  testLocalFiniteElement<C2LBTraits>(&monom7VVFE);

  return 0;

}
catch (Exception e) {

  std::cout << e << std::endl;
  return 1;
}
