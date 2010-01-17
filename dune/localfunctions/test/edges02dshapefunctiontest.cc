// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/grid/common/quadraturerules.hh>
#include <dune/grid/genericgeometry/geometry.hh>

#include <dune/localfunctions/edges02d.hh>

/** \file
    \brief Performs some tests for the monomial shape functions
 */

bool success = true;
double epsilon = 1e-8;

using namespace Dune;

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

// Identity geometry matching reference simplices
template<typename ctype, unsigned dim>
class SimplexGeometry
  : public ReferenceGeometry<ctype, dim>
{
public:
  SimplexGeometry()
    : ReferenceGeometry<ctype, dim>
        (Dune::GeometryType(Dune::GeometryType::simplex,dim))
  {}
};


void testShapeFunctionDerivative(unsigned int orientation, int order)
{

  typedef EdgeS02DLocalFiniteElement<double,double> LocalFiniteElement;
  LocalFiniteElement lFE(orientation);

  typedef LocalFiniteElement::Traits::LocalBasisType LB;

  // ////////////////////////////////////////////////////////////
  //   Check the partial derivatives by comparing them
  //   to finite difference approximations
  // ////////////////////////////////////////////////////////////

  // A set of test points
  // Bad: dependence on dune-grid.  Only for the test points, though
  const QuadratureRule<double,LB::Traits::dimDomain> quad =
    QuadratureRules<double,LB::Traits::dimDomain>::rule(lFE.type(),order);
  const SimplexGeometry<double,LB::Traits::dimDomain> geometry;

  // Loop over all quadrature points
  for (size_t i=0; i<quad.size(); i++) {

    // Get a test point
    const FieldVector<double,LB::Traits::dimDomain>& testPoint = quad[i].position();

    // Get the shape function derivatives there
    std::vector<LB::Traits::JacobianType> jacobians;
    lFE.localBasis().evaluateJacobianGlobal(testPoint, jacobians, geometry);

    // Loop over all shape functions in this set
    for (unsigned int j=0; j<lFE.localBasis().size(); ++j) {

      // Loop over all direction
      for (int k=0; k<LB::Traits::dimDomain; k++) {

        // Compute an approximation to the derivative by finite differences
        FieldVector<double,LB::Traits::dimDomain> upPos   = testPoint;
        FieldVector<double,LB::Traits::dimDomain> downPos = testPoint;

        upPos[k]   += epsilon;
        downPos[k] -= epsilon;

        std::vector<LB::Traits::RangeType> upValues, downValues;

        lFE.localBasis().evaluateFunctionGlobal(upPos,   upValues,   geometry);
        lFE.localBasis().evaluateFunctionGlobal(downPos, downValues, geometry);

        //Loop over all components
        for(int l=0; l < LB::Traits::dimRange; ++l) {

          // The current partial derivative, just for ease of notation
          double derivative = jacobians[j][l][k];

          double finiteDiff = (upValues[j][l] - downValues[j][l]) / (2*epsilon);

          // Check
          if (std::abs(derivative-finiteDiff) > epsilon) {
            std::cerr << "Bug in shape function of order " << order << " for " << lFE.type() << "." << std::endl;
            std::cerr << "Shape function derivative does not agree with FD approximation" << std::endl;
            std::cerr << "Shape function " << j << " component " << l << " at position " << testPoint
                      << ":  derivative in direction " << k << " is " << derivative
                      << ", but " << finiteDiff << " is expected." << std::endl;
            success = false;
          }
        } //Loop over all components
      } // Loop over all direction
    } // Loop over all shape functions in this set
  } // Loop over all quadrature points
}


int main (int argc, char *argv[]) try
{

  testShapeFunctionDerivative(0, 2);

  return success ? 0 : 1;
}
catch (Exception e) {

  std::cout << e << std::endl;
  return 1;
}
