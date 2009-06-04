// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/grid/common/quadraturerules.hh>

#include <dune/finiteelements/monom.hh>

/** \file
    \brief Performs some tests for the monomial shape functions
 */

bool success = true;
double epsilon = 1e-8;

using namespace Dune;

template <int dim, int order>
void testShapeFunctionDerivative(const GeometryType& type)
{

  MonomLocalFiniteElement<double,double,dim,order> shapeFunctionSet(type.basicType());

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

template<int dim, int order>
void testShapeFunctionValue(const GeometryType::BasicType& basicType,
                            const Dune::FieldVector<double, dim> &pos,
                            int comp, double expected)
{
  MonomLocalFiniteElement<double,double,dim,order> shapeFunctionSet(basicType);
  std::vector<Dune::FieldVector<double,1> > out;
  shapeFunctionSet.localBasis().evaluateFunction(pos, out);

  if(std::abs(out[comp][0]-expected) > epsilon) {
    std::cerr << "Bug in shape function of dimension " << dim
              << " and order " << order << " for "
              << GeometryType(basicType,dim) << "." << std::endl;
    std::cerr << "Value of shape function number " << comp << " at position "
              << pos << " is " << out[comp][0] << " but " << expected
              << " was expected." << std::endl;
    success = false;
  }
}

int main (int argc, char *argv[]) try
{
  { // dim=1
    Dune::FieldVector<double, 1> pos;

    pos[0] = 0;
    testShapeFunctionValue<1,2>(GeometryType::cube, pos, 0, 1);
    testShapeFunctionValue<1,2>(GeometryType::cube, pos, 1, 0);
    testShapeFunctionValue<1,2>(GeometryType::cube, pos, 2, 0);

    pos[0] = .5;
    testShapeFunctionValue<1,2>(GeometryType::cube, pos, 0, 1);
    testShapeFunctionValue<1,2>(GeometryType::cube, pos, 1, .5);
    testShapeFunctionValue<1,2>(GeometryType::cube, pos, 2, .25);

    pos[0] = 1;
    testShapeFunctionValue<1,2>(GeometryType::cube, pos, 0, 1);
    testShapeFunctionValue<1,2>(GeometryType::cube, pos, 1, 1);
    testShapeFunctionValue<1,2>(GeometryType::cube, pos, 2, 1);
  }

  { // dim=2
    Dune::FieldVector<double, 2> pos;

    pos[0] = 0; pos[1] = 0;
    testShapeFunctionValue<2,1>(GeometryType::cube, pos, 0, 1);
    testShapeFunctionValue<2,1>(GeometryType::cube, pos, 1, 0);
    testShapeFunctionValue<2,1>(GeometryType::cube, pos, 2, 0);

    pos[0] = .5; pos[1] = .5;
    testShapeFunctionValue<2,1>(GeometryType::cube, pos, 0, 1);
    testShapeFunctionValue<2,1>(GeometryType::cube, pos, 1, .5);
    testShapeFunctionValue<2,1>(GeometryType::cube, pos, 2, .5);

    pos[0] = 1; pos[1] = 1;
    testShapeFunctionValue<2,1>(GeometryType::cube, pos, 0, 1);
    testShapeFunctionValue<2,1>(GeometryType::cube, pos, 1, 1);
    testShapeFunctionValue<2,1>(GeometryType::cube, pos, 2, 1);
  }

  // Test shape functions for the 1d segment
  testShapeFunctionDerivative<1,1>(GeometryType(1));
  testShapeFunctionDerivative<1,2>(GeometryType(1));

  testShapeFunctionDerivative<2,1>(GeometryType(GeometryType::simplex,2));
  //   testShapeFunctionDerivative<2,1>(GeometryType(GeometryType::cube,2));

  //   testShapeFunctionDerivative<3,1>(GeometryType(GeometryType::simplex,3));
  //   testShapeFunctionDerivative<3,1>(GeometryType(GeometryType::cube,3));
  //   testShapeFunctionDerivative<3,1>(GeometryType(GeometryType::pyramid,3));
  //   testShapeFunctionDerivative<3,1>(GeometryType(GeometryType::prism,3));

  //   testShapeFunctionDerivative<4,1>(GeometryType(GeometryType::cube,4));

  return success ? 0 : 1;
}
catch (Exception e) {

  std::cout << e << std::endl;
  return 1;
}
