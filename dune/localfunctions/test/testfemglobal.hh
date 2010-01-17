// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_TESTFEMGLOBAL_HH
#define DUNE_LOCALFUNCTIONS_TESTFEMGLOBAL_HH

#include <iomanip>
#include <iostream>
#include <typeinfo>

#include <dune/grid/common/quadraturerules.hh>
#include <dune/grid/genericgeometry/geometry.hh>
#include <dune/localfunctions/common/virtualinterface.hh>

#include "testfem.hh"

// Identity geometry matching general reference elements
template<typename ctype, unsigned dim>
class ReferenceGeometry
  : public Dune::GenericGeometry::BasicGeometry<
        dim,
        Dune::GenericGeometry::DefaultGeometryTraits<ctype, dim, dim, true> >
{
  typedef Dune::GenericGeometry::BasicGeometry<
      dim,
      Dune::GenericGeometry::DefaultGeometryTraits<ctype, dim, dim, true> >
  Base;
  class RefelemCorners {
    const Dune::GenericReferenceElement<ctype, dim>& refelem;

  public:
    RefelemCorners(const Dune::GeometryType& gt)
      : refelem(Dune::GenericReferenceElements<ctype, dim>::general(gt))
    {}

    const Dune::FieldVector<ctype, dim>&
    operator[](unsigned i) const
    { return refelem.position(i, dim); }
  };

public:
  ReferenceGeometry(const Dune::GeometryType& gt)
    : Base(gt, RefelemCorners(gt))
  {}
};

// This class defines a local finite element function.
// It is determined by a local finite element and
// representing the local basis and a coefficient vector.
// This provides the evaluate method needed by the interpolate()
// method.
template<class FE, typename Geo>
class LocalFEGlobalFunction :
  public Dune::LocalFiniteElementFunctionBase<FE>::type
{
public:
  typedef typename FE::Traits::LocalBasisType::Traits::DomainType DomainType;
  typedef typename FE::Traits::LocalBasisType::Traits::RangeType RangeType;
  typedef typename Dune::Function<const DomainType&, RangeType&> Base;

  typedef typename FE::Traits::LocalBasisType::Traits::RangeFieldType CT;

  LocalFEGlobalFunction(const FE& fe, const Geo& geo_)
    : fe_(fe), geo(geo_)
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
    fe_.localBasis().evaluateFunctionGlobal(x, yy, geo);

    y = 0.0;
    for (std::size_t i=0; i<yy.size(); ++i)
      y.axpy(coeff_[i], yy[i]);
  }

  std::vector<CT> coeff_;

private:
  const FE& fe_;
  const Geo& geo;
};


// Check if localInterpolation is consistens with
// localBasis evaluation.
template<class FE, typename Geo>
bool testLocalInterpolationGlobal(const FE& fe, const Geo& geometry, int n=5)
{
  bool success = true;
  LocalFEGlobalFunction<FE, Geo> f(fe, geometry);

  std::vector<typename LocalFEGlobalFunction<FE, Geo>::CT> coeff;
  for(int i=0; i<n && success; ++i)
  {
    // Set random coefficient vector
    f.setRandom(100);

    // Compute interpolation weights
    fe.localInterpolation().interpolateGlobal(f, coeff, geometry);

    // Check size of weight vector
    if (coeff.size() != fe.localBasis().size())
    {
      std::cout << "Bug in LocalInterpolation for finite element type "
                << typeid(FE).name() << std::endl;
      std::cout << "    Interpolation vector has size " << coeff.size() << std::endl;
      std::cout << "    Basis has size " << fe.localBasis().size() << std::endl;
      std::cout << std::endl;
      success = false;
    }

    // Check if interpolation weights are equal to coefficients
    for(std::size_t j=0; j<coeff.size() && success; ++j)
    {
      if (std::abs(coeff[j]-f.coeff_[j]) > TOL)
      {
        std::cout << std::setprecision(16);
        std::cout << "Bug in LocalInterpolation for finite element type "
                  << typeid(FE).name() << std::endl;
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
template<class FE, typename Geo>
bool testJacobianGlobal(const FE& fe, const Geo& geometry, unsigned order = 2)
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
    fe.localBasis().evaluateJacobianGlobal(testPoint, jacobians, geometry);
    if(jacobians.size() != fe.localBasis().size()) {
      std::cout << "Bug in evaluateJacobianGlobal() for finite element type "
                << typeid(FE).name() << std::endl;
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

        fe.localBasis().evaluateFunctionGlobal(upPos,   upValues,   geometry);
        fe.localBasis().evaluateFunctionGlobal(downPos, downValues, geometry);

        //Loop over all components
        for(int l=0; l < LB::Traits::dimRange; ++l) {

          // The current partial derivative, just for ease of notation
          double derivative = jacobians[j][l][k];

          double finiteDiff = (upValues[j][l] - downValues[j][l])
                              / (2*jacobianTOL);

          // Check
          if (std::abs(derivative-finiteDiff) > jacobianTOL) {
            std::cout << std::setprecision(16);
            std::cout << "Bug in evaluateJacobianGlobal() for finite element "
                      << "type " << typeid(FE).name() << std::endl;
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

// call tests for given finite element
template<class FE>
bool testFEGlobal(const FE& fe, unsigned order = 2)
{
  typedef ReferenceGeometry<
      typename FE::Traits::LocalBasisType::Traits::DomainFieldType,
      FE::Traits::LocalBasisType::Traits::dimDomain> Geo;
  Geo geometry(fe.type());

  std::vector<double> c;

  fe.localInterpolation().interpolateGlobal(Func<FE>(),c, geometry);

  bool success = true;
  success = testLocalInterpolationGlobal(fe, geometry) and success;
  success = testJacobianGlobal(fe, geometry, order) and success;

  return success;
}
#endif // DUNE_LOCALFUNCTIONS_TESTFEMGLOBAL_HH
