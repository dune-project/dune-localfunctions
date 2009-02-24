// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <iostream>

#include <dune/grid/common/quadraturerules.hh>

#include <dune/finiteelements/pk2d.hh>
#include <dune/finiteelements/q12d.hh>

#include <dune/finiteelements/refinedp1.hh>

/** \file
    \brief Performs some tests for the Lagrange shape functions
 */

bool success = true;
double epsilon = 1e-7;

using namespace Dune;

typedef double D;
typedef double R;
enum {k=1};  // Element order

typedef C0LocalBasisTraits<D,2,FieldVector<D,2>,R,1,FieldVector<R,1> > C0Traits;

typedef C1LocalBasisTraits<D,2,FieldVector<D,2>,R,1,FieldVector<R,1>,
    FieldVector<FieldVector<R,2>,1> > C1Traits;


/** \brief Test whether all basis functions in the set sum up to one. */
template <int dim>
void testSumToOne(const C0LocalBasisInterface<C1Traits>* localBasis,
                  const GeometryType& type)
{
  const QuadratureRule<double,dim> quad = QuadratureRules<double,dim>::rule(type,5);

  for (size_t i=0; i<quad.size(); i++) {

    std::vector<FieldVector<double,1> > values;
    localBasis->evaluateFunction(quad[i].position(), values);

    R sum = R(0);
    for (int j=0; j<values.size(); j++)
      sum += values[j];

    if (std::abs(1-sum) > epsilon) {
      std::cerr << "Shapefunction bug: shape functions do not sum up to 1 at " << quad[i].position() << std::endl;
      success = false;
    }

  }

}

/** \brief Test derivative implementations using a finite difference approximation */
template <int dim>
void testShapeFunctionSet(const C1LocalBasisInterface<C1Traits>* localBasis,
                          const GeometryType& type)
{

  // ////////////////////////////////////////////////////////////
  //   Check the partial derivatives by comparing them
  //   to finite difference approximations
  // ////////////////////////////////////////////////////////////

  // A set of test points
  // Bad: dependence on dune-grid.  Only for the test points, though
  const QuadratureRule<double,dim> quad = QuadratureRules<double,dim>::rule(type,5);

  for (size_t i=0; i<quad.size(); i++) {

    // Get a test point
    const FieldVector<double,dim>& testPoint = quad[i].position();

    // Get the shape function derivatives there
    std::vector<FieldVector<FieldVector<double,dim>,1> > jacobians;
    localBasis->evaluateJacobian(testPoint, jacobians);

    // Loop over all shape functions in this set
    for (unsigned int j=0; j<localBasis->size(); ++j) {

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

        localBasis->evaluateFunction(upPos,   upValues);
        localBasis->evaluateFunction(downPos, downValues);
        double finiteDiff = (upValues[j] - downValues[j]) / (2*epsilon);

        // Check
        if (std::abs(derivative-finiteDiff) > epsilon) {
          std::cerr << "Bug in shape function of order " << localBasis->order() << " for " << type << "." << std::endl;
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

  Pk2DLocalFiniteElement<double,double,1> testSetk1;
  testSumToOne<2>(&testSetk1.localBasis(), testSetk1.type());
  testShapeFunctionSet<2>(&testSetk1.localBasis(), testSetk1.type());

  Pk2DLocalFiniteElement<double,double,2> testSetk2;
  testSumToOne<2>(&testSetk2.localBasis(), testSetk2.type());
  testShapeFunctionSet<2>(&testSetk2.localBasis(), testSetk2.type());

  Pk2DLocalFiniteElement<double,double,3> testSetk3;
  testSumToOne<2>(&testSetk3.localBasis(), testSetk3.type());
  testShapeFunctionSet<2>(&testSetk3.localBasis(), testSetk3.type());

  Q12DLocalFiniteElement<double,double> testSetQ1;
  testSumToOne<2>(&testSetQ1.localBasis(), testSetQ1.type());
  testShapeFunctionSet<2>(&testSetQ1.localBasis(), testSetQ1.type());

  RefinedP1LocalFiniteElement<double,double> testSetRefinedP1;
  testSumToOne<2>(&testSetRefinedP1.localBasis(), testSetRefinedP1.type());
  testShapeFunctionSet<2>(&testSetRefinedP1.localBasis(), testSetRefinedP1.type());

  return success ? 0 : 1;
}
catch (Exception e) {

  std::cout << e << std::endl;
  return 1;
}
