// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <config.h>

#include <iostream>

#include <dune/geometry/quadraturerules.hh>

#include <dune/localfunctions/monomial.hh>

/** \file
    \brief Performs some tests for the monomial shape functions
 */

bool success = true;
double epsilon = 1e-8;

using namespace Dune;

template <int dim, int order>
void testShapeFunctionDerivative(const GeometryType& type)
{

  MonomialLocalFiniteElement<double,double,dim,order> shapeFunctionSet(type);
  typedef typename MonomialLocalFiniteElement<double,double,dim,order>::Traits::LocalBasisType::Traits LBTraits;

  // ////////////////////////////////////////////////////////////
  //   Check the partial derivatives by comparing them
  //   to finite difference approximations
  // ////////////////////////////////////////////////////////////

  // A set of test points
  const QuadratureRule<double,dim> quad = QuadratureRules<double,dim>::rule(type,order);

  for (size_t i=0; i<quad.size(); i++) {

    // Get a test point
    const FieldVector<double,dim>& testPoint = quad[i].position();

    // Get the shape function derivatives there
    std::vector<typename LBTraits::JacobianType> jacobians;
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
void testShapeFunctionValue(const GeometryType& gt,
                            const Dune::FieldVector<double, dim> &pos,
                            int comp, double expected)
{
  MonomialLocalFiniteElement<double,double,dim,order> shapeFunctionSet(gt);
  std::vector<Dune::FieldVector<double,1> > out;
  shapeFunctionSet.localBasis().evaluateFunction(pos, out);

  if(std::abs(out[comp][0]-expected) > epsilon) {
    std::cerr << "Bug in shape function of dimension " << dim
              << " and order " << order << " for "
              << gt << "." << std::endl;
    std::cerr << "Value of shape function number " << comp << " at position "
              << pos << " is " << out[comp][0] << " but " << expected
              << " was expected." << std::endl;
    success = false;
  }
}

int main (int argc, char *argv[]) {
  try {
    GeometryType gt;

    {     // dim=1
      Dune::FieldVector<double, 1> pos;
      gt.makeLine();

      pos[0] = 0;
      testShapeFunctionValue<1,2>(gt, pos, 0, 1);
      testShapeFunctionValue<1,2>(gt, pos, 1, 0);
      testShapeFunctionValue<1,2>(gt, pos, 2, 0);

      pos[0] = .5;
      testShapeFunctionValue<1,2>(gt, pos, 0, 1);
      testShapeFunctionValue<1,2>(gt, pos, 1, .5);
      testShapeFunctionValue<1,2>(gt, pos, 2, .25);

      pos[0] = 1;
      testShapeFunctionValue<1,2>(gt, pos, 0, 1);
      testShapeFunctionValue<1,2>(gt, pos, 1, 1);
      testShapeFunctionValue<1,2>(gt, pos, 2, 1);
    }

    {     // dim=2
      Dune::FieldVector<double, 2> pos;
      gt.makeQuadrilateral();

      pos[0] = 0; pos[1] = 0;
      testShapeFunctionValue<2,1>(gt, pos, 0, 1);
      testShapeFunctionValue<2,1>(gt, pos, 1, 0);
      testShapeFunctionValue<2,1>(gt, pos, 2, 0);

      pos[0] = .5; pos[1] = .5;
      testShapeFunctionValue<2,1>(gt, pos, 0, 1);
      testShapeFunctionValue<2,1>(gt, pos, 1, .5);
      testShapeFunctionValue<2,1>(gt, pos, 2, .5);

      pos[0] = 1; pos[1] = 1;
      testShapeFunctionValue<2,1>(gt, pos, 0, 1);
      testShapeFunctionValue<2,1>(gt, pos, 1, 1);
      testShapeFunctionValue<2,1>(gt, pos, 2, 1);
    }

    // Test shape functions for the 1d segment
    gt.makeLine();
    testShapeFunctionDerivative<1,1>(gt);
    testShapeFunctionDerivative<1,2>(gt);

    gt.makeTriangle();
    testShapeFunctionDerivative<2,1>(gt);
    gt.makeQuadrilateral();
    testShapeFunctionDerivative<2,1>(gt);

    gt.makeTetrahedron();
    testShapeFunctionDerivative<3,1>(gt);
    gt.makeHexahedron();
    testShapeFunctionDerivative<3,1>(gt);
    gt.makePyramid();
    testShapeFunctionDerivative<3,1>(gt);
    gt.makePrism();
    testShapeFunctionDerivative<3,1>(gt);

    // gt.makeCube(4);
    // testShapeFunctionDerivative<4,1>(gt);

    return success ? 0 : 1;
  } catch (const Exception &e) {
    std::cout << e << std::endl;
    throw;
  }
}
