// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_TEST_TEST_LOCALFE_HH
#define DUNE_LOCALFUNCTIONS_TEST_TEST_LOCALFE_HH

/** \file \brief Unit tests for LocalFiniteElement objects
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
#include <dune/common/deprecated.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>

double TOL = 1e-9;
// The FD approximation used for checking the Jacobian uses half of the
// precision -- so we have to be a little bit more tolerant here.
double jacobianTOL = 1e-5;  // sqrt(TOL)

template<class FE>
class Func :
  //  public Dune::LocalFiniteElementFunctionBase<FE>::type
  public Dune::LocalFiniteElementFunctionBase<FE>::FunctionBase
  //  public Dune::LocalFiniteElementFunctionBase<FE>::VirtualFunctionBase
{
public:
  typedef typename FE::Traits::LocalBasisType::Traits::DomainType DomainType;
  typedef typename FE::Traits::LocalBasisType::Traits::RangeType RangeType;
  typedef typename Dune::Function<const DomainType&, RangeType&> Base;

  void evaluate (const DomainType& x, RangeType& y) const
  {
    y = 0;
    DomainType c(0.5);

    c -= x;
    y[0] = exp(-3.0*c.two_norm2());
  }
};

// This class defines a local finite element function.
// It is determined by a local finite element and
// representing the local basis and a coefficient vector.
// This provides the evaluate method needed by the interpolate()
// method.
template<class FE>
class LocalFEFunction :
  //  public Dune::LocalFiniteElementFunctionBase<FE>::type
  public Dune::LocalFiniteElementFunctionBase<FE>::FunctionBase
  //  public Dune::LocalFiniteElementFunctionBase<FE>::VirtualFunctionBase
{
public:
  typedef typename FE::Traits::LocalBasisType::Traits::DomainType DomainType;
  typedef typename FE::Traits::LocalBasisType::Traits::RangeType RangeType;
  typedef typename Dune::Function<const DomainType&, RangeType&> Base;

  typedef typename FE::Traits::LocalBasisType::Traits::RangeFieldType CT;

  LocalFEFunction(const FE& fe) :
    fe_(fe)
  {
    resetCoefficients();
  }

  void resetCoefficients()
  {
    coeff_.resize(fe_.localBasis().size());
    for(std::size_t i=0; i<coeff_.size(); ++i)
      coeff_[i] = 0;
  }

  void setRandom(double max)
  {
    coeff_.resize(fe_.localBasis().size());
    for(std::size_t i=0; i<coeff_.size(); ++i)
      coeff_[i] = ((1.0*std::rand()) / RAND_MAX - 0.5)*2.0*max;
  }


  void evaluate (const DomainType& x, RangeType& y) const
  {
    std::vector<RangeType> yy;
    fe_.localBasis().evaluateFunction(x, yy);

    y = 0.0;
    for (std::size_t i=0; i<yy.size(); ++i)
      y.axpy(coeff_[i], yy[i]);
  }

  std::vector<CT> coeff_;

private:
  const FE& fe_;
};


// Check if localInterpolation is consistens with
// localBasis evaluation.
template<class FE>
bool testLocalInterpolation(const FE& fe, int n=5)
{
  bool success = true;
  LocalFEFunction<FE> f(fe);

  std::vector<typename LocalFEFunction<FE>::CT> coeff;
  for(int i=0; i<n && success; ++i)
  {
    // Set random coefficient vector
    f.setRandom(100);

    // Compute interpolation weights
    fe.localInterpolation().interpolate(f, coeff);

    // Check size of weight vector
    if (coeff.size() != fe.localBasis().size())
    {
      std::cout << "Bug in LocalInterpolation for finite element type "
                << Dune::className(fe) << std::endl;
      std::cout << "    Interpolation vector has size " << coeff.size() << std::endl;
      std::cout << "    Basis has size " << fe.localBasis().size() << std::endl;
      std::cout << std::endl;
      success = false;
    }

    // Check if interpolation weights are equal to coefficients
    for(std::size_t j=0; j<coeff.size() && success; ++j)
    {
      if ( std::abs(coeff[j]-f.coeff_[j]) >
           TOL*((std::abs(f.coeff_[j])>1) ? std::abs(f.coeff_[j]) : 1.) )
      {
        std::cout << std::setprecision(16);
        std::cout << "Bug in LocalInterpolation for finite element type "
                  << Dune::className(fe) << std::endl;
        std::cout << "    Interpolation weight " << j
                  << " differs by " << std::abs(coeff[j]-f.coeff_[j])
                  << " from coefficient of linear combination." << std::endl;
        std::cout << std::endl;
        success = false;
      }
    }
  }
  return success;
}


// check whether Jacobian agrees with FD approximation
template<class FE>
bool testJacobian(const FE& fe, unsigned order = 2)
{
  typedef typename FE::Traits::LocalBasisType LB;

  bool success = true;

  // ////////////////////////////////////////////////////////////
  //   Check the partial derivatives by comparing them
  //   to finite difference approximations
  // ////////////////////////////////////////////////////////////

  // A set of test points
  const Dune::QuadratureRule<double,LB::Traits::dimDomain> quad =
    Dune::QuadratureRules<double,LB::Traits::dimDomain>::rule(fe.type(),order);

  // Loop over all quadrature points
  for (size_t i=0; i<quad.size(); i++) {

    // Get a test point
    const Dune::FieldVector<double,LB::Traits::dimDomain>& testPoint =
      quad[i].position();

    // Get the shape function derivatives there
    std::vector<typename LB::Traits::JacobianType> jacobians;
    fe.localBasis().evaluateJacobian(testPoint, jacobians);
    if(jacobians.size() != fe.localBasis().size()) {
      std::cout << "Bug in evaluateJacobianGlobal() for finite element type "
                << Dune::className(fe) << std::endl;
      std::cout << "    Jacobian vector has size " << jacobians.size()
                << std::endl;
      std::cout << "    Basis has size " << fe.localBasis().size()
                << std::endl;
      std::cout << std::endl;
      return false;
    }

    // Loop over all shape functions in this set
    for (unsigned int j=0; j<fe.localBasis().size(); ++j) {
      // Loop over all directions
      for (int k=0; k<LB::Traits::dimDomain; k++) {

        // Compute an approximation to the derivative by finite differences
        Dune::FieldVector<double,LB::Traits::dimDomain> upPos   = testPoint;
        Dune::FieldVector<double,LB::Traits::dimDomain> downPos = testPoint;

        upPos[k]   += jacobianTOL;
        downPos[k] -= jacobianTOL;

        std::vector<typename LB::Traits::RangeType> upValues, downValues;

        fe.localBasis().evaluateFunction(upPos,   upValues);
        fe.localBasis().evaluateFunction(downPos, downValues);

        //Loop over all components
        for(int l=0; l < LB::Traits::dimRange; ++l) {

          // The current partial derivative, just for ease of notation
          double derivative = jacobians[j][l][k];

          double finiteDiff = (upValues[j][l] - downValues[j][l])
                              / (2*jacobianTOL);

          // Check
          if ( std::abs(derivative-finiteDiff) >
               TOL/jacobianTOL*((std::abs(finiteDiff)>1) ? std::abs(finiteDiff) : 1.) )
          {
            std::cout << std::setprecision(16);
            std::cout << "Bug in evaluateJacobian() for finite element type "
                      << Dune::className(fe) << std::endl;
            std::cout << "    Shape function derivative does not agree with "
                      << "FD approximation" << std::endl;
            std::cout << "    Shape function " << j << " component " << l
                      << " at position " << testPoint << ": derivative in "
                      << "direction " << k << " is " << derivative << ", but "
                      << finiteDiff << " is expected." << std::endl;
            std::cout << std::endl;
            success = false;
          }
        } //Loop over all components
      } // Loop over all directions
    } // Loop over all shape functions in this set
  } // Loop over all quadrature points

  return success;
}

/** \brief Helper class to test the 'evaluate' method
 *
 * It implements a static loop over the available diff orders
 */
template<int diffOrder>
struct TestEvaluate
{
  template <class FE>
  static bool test(const FE& fe,
                   double eps, double delta, std::size_t order = 2)
  {
    std::cout << "No test for differentiability order " << diffOrder << std::endl;
    return TestEvaluate<diffOrder-1>::test(fe, eps, delta, order);
  }
};

/** \brief Specialization to test the 'evaluate' method for zero-order partial derivatives, i.e., values */
template<>
struct TestEvaluate<0>
{
  template <class FE>
  static bool test(const FE& fe,
                   double eps,
                   double delta,
                   std::size_t order = 2)
  {
    typedef typename FE::Traits::LocalBasisType::Traits::RangeType RangeType;
    constexpr auto dimDomain = FE::Traits::LocalBasisType::Traits::dimDomain;

    bool success = true;

    //////////////////////////////////////////////////////////////
    //   Check the partial derivatives by comparing them
    //   to finite difference approximations
    //////////////////////////////////////////////////////////////

    // A set of test points
    const auto& quad = Dune::QuadratureRules<double, dimDomain>::rule(fe.type(),
                                                                      order);

    // Loop over all quadrature points
    for (size_t i = 0; i < quad.size(); i++)
    {
      // Get a test point
      const Dune::FieldVector<double, dimDomain>& testPoint = quad[i].position();

      // Get the shape function values there using the 'partial' method
      std::vector<RangeType> partialValues;
      std::array<unsigned int, dimDomain> multiIndex;
      std::fill(multiIndex.begin(), multiIndex.end(), 0);

      fe.localBasis().partial(multiIndex, testPoint, partialValues);

      if (partialValues.size() != fe.localBasis().size())
      {
        std::cout << "Bug in partial() for finite element type "
                  << Dune::className(fe) << std::endl;
        std::cout << "    values vector has size "
                  << partialValues.size() << std::endl;
        std::cout << "    Basis has size " << fe.localBasis().size()
                  << std::endl;
        std::cout << std::endl;
        return false;
      }

      // Get reference values
      std::vector<RangeType> referenceValues;
      fe.localBasis().evaluateFunction(testPoint, referenceValues);

      // Loop over all shape functions in this set
      for (unsigned int j = 0; j < fe.localBasis().size(); ++j)
      {
        // Loop over all components
        for (int l = 0; l < FE::Traits::LocalBasisType::Traits::dimRange; ++l)
        {
          // Check the 'partial' method
          if (std::abs(partialValues[j][l] - referenceValues[j][l])
              > TOL / jacobianTOL
                * ((std::abs(referenceValues[j][l]) > 1) ? std::abs(referenceValues[j][l]) : 1.))
          {
            std::cout << std::setprecision(16);
            std::cout << "Bug in partial() for finite element type "
                      << Dune::className(fe) << std::endl;
            std::cout << "    Shape function value does not agree with "
                      << "output of method evaluateFunction." << std::endl;
            std::cout << "    Shape function " << j << " component " << l
                      << " at position " << testPoint << ": value is " << partialValues[j][l]
                      << ", but " << referenceValues[j][l] << " is expected." << std::endl;
            std::cout << std::endl;
            success = false;
          }

        } //Loop over all components
      } // Loop over all shape functions in this set
    } // Loop over all quadrature points

    return success;
  }
};

/** \brief Specialization to test the 'evaluate' method for first-order partial derivatives */
template<>
struct TestEvaluate<1>
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

      // Loop over all directions
      for (int k = 0; k < LB::Traits::dimDomain; k++)
      {
        std::array<int, 1> direction = {{k}};

        // Get the shape function derivatives there using the 'evaluate' method
        std::vector<typename LB::Traits::RangeType> firstDerivatives;
        DUNE_NO_DEPRECATED_BEGIN
        fe.localBasis().template evaluate<1>(direction, testPoint, firstDerivatives);
        DUNE_NO_DEPRECATED_END
        if (firstDerivatives.size() != fe.localBasis().size())
        {
          std::cout << "Bug in evaluate() for finite element type "
                    << Dune::className(fe) << std::endl;
          std::cout << "    firstDerivatives vector has size "
                    << firstDerivatives.size() << std::endl;
          std::cout << "    Basis has size " << fe.localBasis().size()
                    << std::endl;
          std::cout << std::endl;
          return false;
        }

        // Get the shape function derivatives there using the 'partial' method
        std::vector<typename LB::Traits::RangeType> firstPartialDerivatives;
        std::array<unsigned int, LB::Traits::dimDomain> multiIndex;
        std::fill(multiIndex.begin(), multiIndex.end(), 0);
        multiIndex[k]++;
        fe.localBasis().partial(multiIndex, testPoint, firstPartialDerivatives);
        if (firstPartialDerivatives.size() != fe.localBasis().size())
        {
          std::cout << "Bug in evaluate() for finite element type "
                    << Dune::className(fe) << std::endl;
          std::cout << "    firstDerivatives vector has size "
                    << firstDerivatives.size() << std::endl;
          std::cout << "    Basis has size " << fe.localBasis().size()
                    << std::endl;
          std::cout << std::endl;
          return false;
        }

        // Loop over all shape functions in this set
        for (unsigned int j = 0; j < fe.localBasis().size(); ++j)
        {
          // Compute an approximation to the derivative by finite differences
          Dune::FieldVector<double, LB::Traits::dimDomain> upPos = testPoint;
          Dune::FieldVector<double, LB::Traits::dimDomain> downPos = testPoint;

          upPos[k] += jacobianTOL;
          downPos[k] -= jacobianTOL;

          std::vector<typename LB::Traits::RangeType> upValues, downValues;

          fe.localBasis().evaluateFunction(upPos, upValues);
          fe.localBasis().evaluateFunction(downPos, downValues);

          // Loop over all components
          for (int l = 0; l < LB::Traits::dimRange; ++l)
          {
            // The current partial derivative, just for ease of notation
            RangeField derivative = firstDerivatives[j][l];

            RangeField finiteDiff = (upValues[j][l] - downValues[j][l])
                              / (2 * jacobianTOL);

            // Check the 'evaluate' method
            if (std::abs(derivative - finiteDiff)
                > TOL / jacobianTOL
                  * ((std::abs(finiteDiff) > 1) ? std::abs(finiteDiff) : 1.))
            {
              std::cout << std::setprecision(16);
              std::cout << "Bug in evaluate<1>() for finite element type "
                        << Dune::className(fe) << std::endl;
              std::cout << "    Shape function derivative does not agree with "
                        << "FD approximation" << std::endl;
              std::cout << "    Shape function " << j << " component " << l
                        << " at position " << testPoint << ": derivative in "
                        << "direction " << k << " is " << derivative << ", but "
                        << finiteDiff << " is expected." << std::endl;
              std::cout << std::endl;
              success = false;
            }

            // Check the 'partial' method
            RangeField partialDerivative = firstPartialDerivatives[j][l];
            if (std::abs(partialDerivative - finiteDiff)
                > TOL / jacobianTOL
                  * ((std::abs(finiteDiff) > 1) ? std::abs(finiteDiff) : 1.))
            {
              std::cout << std::setprecision(16);
              std::cout << "Bug in partial() for finite element type "
                        << Dune::className(fe) << std::endl;
              std::cout << "    Shape function derivative does not agree with "
                        << "FD approximation" << std::endl;
              std::cout << "    Shape function " << j << " component " << l
                        << " at position " << testPoint << ": derivative in "
                        << "direction " << k << " is " << derivative << ", but "
                        << finiteDiff << " is expected." << std::endl;
              std::cout << std::endl;
              success = false;
            }

          } // Loop over all directions
        } //Loop over all components
      } // Loop over all shape functions in this set
    } // Loop over all quadrature points

    // Recursively call the zero-order test
    return success and TestEvaluate<0>::test(fe, eps, delta, order);
  }
};

/** \brief Specialization to test second-order partial derivatives */
template<>
struct TestEvaluate<2>
{
  template <class FE>
  static bool test(const FE& fe,
                   double eps,
                   double delta,
                   std::size_t order = 2)
  {
    typedef typename FE::Traits::LocalBasisType LocalBasis;
    typedef typename LocalBasis::Traits::DomainFieldType DF;
    typedef typename LocalBasis::Traits::DomainType Domain;
    static const int dimDomain = LocalBasis::Traits::dimDomain;

    static const std::size_t dimR = LocalBasis::Traits::dimRange;
    typedef typename LocalBasis::Traits::RangeType Range;
    typedef typename LocalBasis::Traits::RangeFieldType RangeField;

    bool success = true;

    //////////////////////////////////////////////////////////////
    //   Check the partial derivatives by comparing them
    //   to finite difference approximations
    //////////////////////////////////////////////////////////////

    // A set of test points
    const Dune::QuadratureRule<DF, dimDomain> quad
       = Dune::QuadratureRules<DF,dimDomain>::rule(fe.type(), order);

    // Loop over all quadrature points
    for (std::size_t i = 0; i < quad.size(); i++)
    {
      // Get a test point
      const Domain& testPoint = quad[i].position();

      // For testing the 'evaluate' method
      std::array<std::vector<Dune::FieldMatrix<RangeField, dimDomain, dimDomain> >, dimR> hessians;
      for (size_t k = 0; k < dimR; k++)
        hessians[k].resize(fe.size());

      // For testing the 'partial' method
      std::array<std::vector<Dune::FieldMatrix<RangeField, dimDomain, dimDomain> >, dimR> partialHessians;
      for (size_t k = 0; k < dimR; k++)
        partialHessians[k].resize(fe.size());

      //loop over all local directions
      for (int dir0 = 0; dir0 < dimDomain; dir0++)
      {
        for (int dir1 = 0; dir1 < dimDomain; dir1++)
        {
          std::array<int, 2> directions = {{ dir0, dir1 }};

          // Get the shape function derivatives there using the 'evaluate' method
          std::vector<Range> secondDerivative;
          DUNE_NO_DEPRECATED_BEGIN
          fe.localBasis().template evaluate<2>(directions, testPoint, secondDerivative);
          DUNE_NO_DEPRECATED_END
          if (secondDerivative.size() != fe.localBasis().size())
          {
            std::cout << "Bug in evaluate<2>() for finite element type "
                      << Dune::className<FE>() << ":" << std::endl;
            std::cout << "    return vector has size " << secondDerivative.size()
                      << std::endl;
            std::cout << "    Basis has size " << fe.localBasis().size()
                      << std::endl;
            std::cout << std::endl;
            return false;
          }

          //combine to Hesse matrices
          for (size_t k = 0; k < dimR; k++)
            for (std::size_t j = 0; j < fe.localBasis().size(); ++j)
              hessians[k][j][dir0][dir1] = secondDerivative[j][k];

          // Get the shape function derivatives there using the 'partial' method
          std::vector<Range> secondPartialDerivative;
          std::array<unsigned int,dimDomain> multiIndex;
          std::fill(multiIndex.begin(), multiIndex.end(), 0);
          multiIndex[dir0]++;
          multiIndex[dir1]++;
          fe.localBasis().partial(multiIndex, testPoint, secondPartialDerivative);

          if (secondPartialDerivative.size() != fe.localBasis().size())
          {
            std::cout << "Bug in partial() for finite element type "
                      << Dune::className<FE>() << ":" << std::endl;
            std::cout << "    return vector has size " << secondPartialDerivative.size()
                      << std::endl;
            std::cout << "    Basis has size " << fe.localBasis().size()
                      << std::endl;
            std::cout << std::endl;
            return false;
          }

          //combine to Hesse matrices
          for (size_t k = 0; k < dimR; k++)
            for (std::size_t j = 0; j < fe.localBasis().size(); ++j)
              partialHessians[k][j][dir0][dir1] = secondPartialDerivative[j][k];

        }
      }  //loop over all directions

      // Loop over all shape functions in this set
      for (std::size_t j = 0; j < fe.localBasis().size(); ++j)
      {
        // Loop over all local directions
        for (std::size_t dir0 = 0; dir0 < dimDomain; ++dir0)
        {
          for (unsigned int dir1 = 0; dir1 < dimDomain; dir1++)
          {
            // Compute an approximation to the derivative by finite differences
            std::array<Domain,4> neighbourPos;
            std::fill(neighbourPos.begin(), neighbourPos.end(), testPoint);

            neighbourPos[0][dir0] += delta;
            neighbourPos[0][dir1] += delta;
            neighbourPos[1][dir0] -= delta;
            neighbourPos[1][dir1] += delta;
            neighbourPos[2][dir0] += delta;
            neighbourPos[2][dir1] -= delta;
            neighbourPos[3][dir0] -= delta;
            neighbourPos[3][dir1] -= delta;

            std::array<std::vector<Range>, 4> neighbourValues;
            for (int k = 0; k < 4; k++)
              fe.localBasis().evaluateFunction(neighbourPos[k],
                                               neighbourValues[k]);

            //Loop over all components
            for (std::size_t k = 0; k < dimR; ++k)
            {
              // The current partial derivative, just for ease of notation, evaluated by the 'evaluate' method
              RangeField derivative = hessians[k][j][dir0][dir1];

              RangeField finiteDiff = (neighbourValues[0][j][k]
                  - neighbourValues[1][j][k] - neighbourValues[2][j][k]
                  + neighbourValues[3][j][k]) / (4 * delta * delta);

              // Check
              if (std::abs(derivative - finiteDiff)
                  > eps / delta * (std::max(std::abs(finiteDiff), 1.0)))
              {
                std::cout << std::setprecision(16);
                std::cout << "Bug in evaluate<2>() for finite element type "
                          << Dune::className<FE>() << ":" << std::endl;
                std::cout << "    Second shape function derivative does not agree with "
                          << "FD approximation" << std::endl;
                std::cout << "    Shape function " << j << " component " << k
                          << " at position " << testPoint << ": derivative in "
                          << "local direction (" << dir0 << ", " << dir1 << ") is "
                          << derivative << ", but " << finiteDiff
                          << " is expected." << std::endl;
                std::cout << std::endl;
                success = false;
              }

              // The current partial derivative, just for ease of notation, evaluated by the 'partial' method
              RangeField partialDerivative = partialHessians[k][j][dir0][dir1];

              // Check
              if (std::abs(partialDerivative - finiteDiff)
                  > eps / delta * (std::max(std::abs(finiteDiff), 1.0)))
              {
                std::cout << std::setprecision(16);
                std::cout << "Bug in partial() for finite element type "
                          << Dune::className<FE>() << ":" << std::endl;
                std::cout << "    Second shape function derivative does not agree with "
                          << "FD approximation" << std::endl;
                std::cout << "    Shape function " << j << " component " << k
                          << " at position " << testPoint << ": derivative in "
                          << "local direction (" << dir0 << ", " << dir1 << ") is "
                          << partialDerivative << ", but " << finiteDiff
                          << " is expected." << std::endl;
                std::cout << std::endl;
                success = false;
              }

            } //Loop over all components
          }
        } // Loop over all local directions
      } // Loop over all shape functions in this set
    } // Loop over all quadrature points

    // Recursively call the first-order test
    return success and TestEvaluate<1>::test(fe, eps, delta, order);
  }

};


template<>
struct TestEvaluate<-1>
{
  template <class FE>
  static bool test(const FE& fe,
                   double eps,
                   double delta,
                   std::size_t order = 2)
  {
    return true;
  }
};

// Flags for disabling parts of testFE
enum {
  DisableNone = 0,
  DisableLocalInterpolation = 1,
  DisableVirtualInterface = 2,
  DisableJacobian = 4,
  DisableEvaluate = 8
};

// call tests for given finite element
template<class FE>
bool testFE(const FE& fe, char disabledTests = DisableNone, unsigned order = 2)
{
  std::vector<double> c;


  bool success = true;

  if (FE::Traits::LocalBasisType::Traits::dimDomain != fe.type().dim())
  {
    std::cout << "Bug in type() for finite element type "
              << Dune::className(fe) << std::endl;
    std::cout << "    Coordinate dimension is " << FE::Traits::LocalBasisType::Traits::dimDomain << std::endl;
    std::cout << "    but GeometryType is " << fe.type() << " with dimension " << fe.type().dim() << std::endl;
    success = false;
  }

  if (fe.size() != fe.localBasis().size())
  {
    std::cout << "Bug in finite element type "
              << Dune::className(fe) << std::endl;
    std::cout << "    Size reported by LocalFiniteElement is " << fe.size() << std::endl;
    std::cout << "    but size reported by LocalBasis is " << fe.localBasis().size() << std::endl;
    success = false;
  }

  // Make sure evaluateFunction returns the correct number of values
  std::vector<typename FE::Traits::LocalBasisType::Traits::RangeType> values;
  fe.localBasis().evaluateFunction(Dune::ReferenceElements<double,FE::Traits::LocalBasisType::Traits::dimDomain>::general(fe.type()).position(0,0), values);

  if (values.size() != fe.size())
  {
    std::cout << "Bug in finite element type "
              << Dune::className(fe) << std::endl;
    std::cout << "    LocalFiniteElement.size() returns " << fe.size() << "," << std::endl;
    std::cout << "    but LocalBasis::evaluateFunction returns " << values.size() << " values!" << std::endl;
    success = false;
  }

  if (fe.size() != fe.localCoefficients().size())
  {
    std::cout << "Bug in finite element type "
              << Dune::className(fe) << std::endl;
    std::cout << "    Size reported by LocalFiniteElement is " << fe.size() << std::endl;
    std::cout << "    but size reported by LocalCoefficients is " << fe.localCoefficients().size() << std::endl;
    success = false;
  }

  if (not (disabledTests & DisableLocalInterpolation))
  {
    fe.localInterpolation().interpolate(Func<FE>(),c);
    success = testLocalInterpolation<FE>(fe) and success;
  }
  if (not (disabledTests & DisableJacobian))
  {
    success = testJacobian<FE>(fe, order) and success;
  }
  else
  {
    // make sure diffOrder is 0
    success = (FE::Traits::LocalBasisType::Traits::diffOrder == 0) and success;
  }

  if (not (disabledTests & DisableEvaluate))
  {
    success = TestEvaluate<FE::Traits::LocalBasisType::Traits::diffOrder>::test(fe, TOL, jacobianTOL, order) and success;
  }

  if (not (disabledTests & DisableVirtualInterface))
  {
    typedef typename FE::Traits::LocalBasisType::Traits LBTraits;
    typedef typename Dune::FixedOrderLocalBasisTraits<LBTraits,0>::Traits C0LBTraits;
    typedef typename Dune::LocalFiniteElementVirtualInterface<C0LBTraits> VirtualFEInterface;
    typedef typename Dune::LocalFiniteElementVirtualImp<FE> VirtualFEImp;

    const VirtualFEImp virtualFE(fe);
    if (not (disabledTests & DisableLocalInterpolation))
      success = testLocalInterpolation<VirtualFEInterface>(virtualFE) and success;
    if (not (disabledTests & DisableJacobian))
    {
      success = testJacobian<VirtualFEInterface>(virtualFE) and success;
    }
    else
    {
      // make sure diffOrder is 0
      success = (VirtualFEInterface::Traits::LocalBasisType::Traits::diffOrder == 0) and success;
    }
  }

  return success;
}

#define TEST_FE(A) { bool b = testFE(A); std::cout << "testFE(" #A ") " << (b?"succeeded\n":"failed\n"); success &= b; }
#define TEST_FE2(A,B) { bool b = testFE(A, B); if (!b) std::cerr << "testFE(" #A ", " #B ") " << (b?"succeeded\n":"failed\n"); success &= b; }

#endif // DUNE_LOCALFUNCTIONS_TEST_TEST_LOCALFE_HH
