// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_TEST_TEST_LOCALFE_HH
#define DUNE_LOCALFUNCTIONS_TEST_TEST_LOCALFE_HH

/** \file \brief Unit tests for LocalFiniteElement objects with pyramid geometry
 *
 * \note This header is not part of the official Dune API and might be subject
 *  to change.  You can use this header to test external finite element
 *  implementations, but be warned that your tests might break with future
 *  Dune versions.
 */

#include <iomanip>
#include <iostream>
#include <typeinfo>

#include <dune/common/classname.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/common/partial.hh>
#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>

double TOL = 1e-9;
// The FD approximation used for checking the Jacobian uses half of the
// precision -- so we have to be a little bit more tolerant here.
double jacobianTOL = 1e-5;  // sqrt(TOL)

template <class T,
  class = std::enable_if_t< std::is_arithmetic<T>::value > >
T sqr(T const& x)
{
  return x*x;
}

template <class T,
  class = std::enable_if_t< std::is_arithmetic<T>::value > >
T distance(T const& x, T const& y)
{
  return std::abs(x - y);
}

template <class T, int n,
  class = std::enable_if_t< std::is_arithmetic<T>::value > >
T distance(Dune::FieldVector<T, n> const& x, Dune::FieldVector<T, n> const& y)
{
  T result = 0;
  for (int i = 0; i < n; ++i)
    result += sqr(x[i] - y[i]);
  return std::sqrt(result);
}


/**
 * \brief Specialization to test the 'evaluate' method for first-order partial derivatives
 * Compare the partial() derivative calculation agains evaluateJacobian() method for
 * consistency.
 **/
struct TestEvaluate
{
  template <class FE>
  static bool test(const FE& fe,
                   double eps,
                   double delta,
                   std::size_t order = 2)
  {
    typedef typename FE::Traits::LocalBasisType LB;
    typedef typename LB::Traits::RangeFieldType RangeField;

    bool success = true;

    //////////////////////////////////////////////////////////////
    //   Check the partial derivatives by comparing them
    //   to finite difference approximations
    //////////////////////////////////////////////////////////////

    // A set of test points
    const Dune::QuadratureRule<double, LB::Traits::dimDomain> quad =
      Dune::QuadratureRules<double, LB::Traits::dimDomain>::rule(fe.type(),
                                                                 order);

    // Loop over all quadrature points
    for (size_t i = 0; i < quad.size(); i++)
    {
      // Get a test point
      const Dune::FieldVector<double, LB::Traits::dimDomain>& testPoint = quad[i].position();


      std::vector<typename LB::Traits::JacobianType> gradients;
      fe.localBasis().evaluateJacobian(testPoint, gradients);


      // Loop over all directions
      for (int k = 0; k < LB::Traits::dimDomain; k++)
      {
        std::array<int, 1> direction = {{k}};


        // Get the shape function derivatives there using the 'partial' method
        std::vector<typename LB::Traits::RangeType> partialDerivatives;

        std::array<unsigned int, LB::Traits::dimDomain> multiIndex;
        Dune::Impl::directions2order(direction, multiIndex);

        fe.localBasis().partial(multiIndex, testPoint, partialDerivatives);
        if (partialDerivatives.size() != gradients.size())
        {
          std::cout << "Bug in evaluate() for finite element type "
                    << Dune::className(fe) << std::endl;
          std::cout << "    partialDerivatives vector has size "
                    << partialDerivatives.size() << std::endl;
          std::cout << "    Basis has size " << fe.localBasis().size()
                    << std::endl;
          std::cout << std::endl;
          return false;
        }

        // Loop over all shape functions in this set
        for (unsigned int j = 0; j < fe.localBasis().size(); ++j)
        {
          // Loop over all components
          for (int l = 0; l < LB::Traits::dimRange; ++l)
          {
            // The current partial derivative, just for ease of notation
            RangeField derivative = partialDerivatives[j][l];

            // Check the 'evaluate' method
            if (distance(derivative, gradients[j][l][k])
                > TOL * std::abs(gradients[j][l][k]))
            {
              std::cout << std::setprecision(16);
              std::cout << "Bug in partial() for finite element type "
                        << Dune::className(fe) << std::endl;
              std::cout << "    Shape function derivative does not agree with "
                        << "FD approximation" << std::endl;
              std::cout << "    Shape function " << j << " component " << l
                        << " at position " << testPoint << ": derivative in "
                        << "direction " << k << " is " << derivative << ", but "
                        << gradients[j][l][k] << " is expected." << std::endl;
              std::cout << std::endl;
              success = false;
            }

          } // Loop over all directions
        } //Loop over all components
      } // Loop over all shape functions in this set
    } // Loop over all quadrature points

    return success;
  }
};


// call tests for given finite element
template<class FE>
bool testFE(const FE& fe)
{
  std::vector<double> c;

  bool success = true;
  success = TestEvaluate::test(fe, TOL, jacobianTOL) and success;

  return success;
}

#endif // DUNE_LOCALFUNCTIONS_TEST_TEST_LOCALFE_HH
