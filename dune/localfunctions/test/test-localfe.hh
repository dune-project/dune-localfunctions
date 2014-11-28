// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_TEST_TEST_LOCALFE_HH
#define DUNE_LOCALFUNCTIONS_TEST_TEST_LOCALFE_HH

#include <iomanip>
#include <iostream>
#include <typeinfo>

#include <dune/common/classname.hh>
#include <dune/geometry/quadraturerules.hh>

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


// Flags for disabling parts of testFE
enum {
  DisableNone = 0,
  DisableLocalInterpolation = 1,
  DisableVirtualInterface = 2,
  DisableJacobian = 3
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
#define TEST_FE2(A,B) { bool b = testFE(A, B); if (!b) std::cerr << "testFE(" #A ", " #B ") " << (b?"succeded\n":"failed\n"); success &= b; }

#endif // DUNE_LOCALFUNCTIONS_TEST_TEST_LOCALFE_HH
