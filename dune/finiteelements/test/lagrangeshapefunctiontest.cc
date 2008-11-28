// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/grid/common/quadraturerules.hh>

#include <dune/finiteelements/pk2d.hh>

/** \file
    \brief Performs some tests for the Lagrange shape functions
 */

bool success = true;
double epsilon = 1e-8;

using namespace Dune;

template <int dim>
void testShapeFunction(const GeometryType& type, int order)
{

  Pk2DLocalFiniteElement<double,double,1> shapeFunctionSet;

  // ////////////////////////////////////////////////////////////
  //   Check the partial derivatives by comparing them
  //   to finite difference approximations
  // ////////////////////////////////////////////////////////////

  // A set of test points
  // Bad: dependence on dune-grid.  Only for the test points, though
  const QuadratureRule<double,dim> quad = QuadratureRules<double,dim>::rule(type,order);

  for (size_t i=0; i<quad.size(); i++) {

    // Get a test point
    const FieldVector<double,dim>& testPoint = quad[i].position();

    // Get the shape function derivatives there
    std::vector<FieldVector<FieldVector<double,dim>,1> > jacobians;
    shapeFunctionSet.localBasis().evaluateJacobian(testPoint, jacobians);

    // Loop over all shape functions in this set
    for (unsigned int j=0; j<shapeFunctionSet.localBasis().size(); ++j) {

      // Loop over all partial derivatives
      for (int k=0; k<dim; k++) {

        // The current partial derivative, just for ease of notation
        double derivative = jacobians[j][0][k];

        // Compute an approximation to the derivative by finite differences
        FieldVector<double,dim> upPos   = testPoint;
        FieldVector<double,dim> downPos = testPoint;

        upPos[k]   += epsilon;
        downPos[k] -= epsilon;

        std::vector<FieldVector<double,1> > upValues, downValues;

        shapeFunctionSet.localBasis().evaluateFunction(upPos,   upValues);
        shapeFunctionSet.localBasis().evaluateFunction(downPos, downValues);
        double finiteDiff = (upValues[j] - downValues[j]) / (2*epsilon);

        // Check
        if (std::abs(derivative-finiteDiff) > epsilon) {
          std::cerr << "Bug in shape function of order " << order << " for " << type << "." << std::endl;
          std::cerr << "Shape function derivative does not agree with FD approximation" << std::endl;
          std::cerr << "Shape function " << j << " at position " << testPoint
                    << ":  derivative in direction " << k << " is " << derivative
                    << ", but " << finiteDiff << " is expected." << std::endl;
          success = false;
        }

      }

    }

  }

}

int main (int argc, char *argv[]) try
{
  for (int i=0; i<3; i++) {

#if 0
    // Test shape functions for the 1d segment
    testShapeFunction<1>(GeometryType(1),i);
#endif
    testShapeFunction<2>(GeometryType(GeometryType::simplex,2),i);
#if 0
    testShapeFunction<2>(GeometryType(GeometryType::cube,2),i);

    testShapeFunction<3>(GeometryType(GeometryType::simplex,3),i);
    testShapeFunction<3>(GeometryType(GeometryType::cube,3),i);
    testShapeFunction<3>(GeometryType(GeometryType::pyramid,3),i);
    testShapeFunction<3>(GeometryType(GeometryType::prism,3),i);

    testShapeFunction<4>(GeometryType(GeometryType::cube,4),i);
#endif
  }

  return success ? 0 : 1;
}
catch (Exception e) {

  std::cout << e << std::endl;
  return 1;
}
