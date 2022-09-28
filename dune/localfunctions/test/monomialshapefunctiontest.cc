// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#include <config.h>

#include <iostream>

#include <dune/geometry/quadraturerules.hh>

#include <dune/localfunctions/monomial.hh>
#include <dune/localfunctions/test/test-localfe.hh>

/** \file
    \brief Performs some tests for the monomial shape functions
 */

double epsilon = 1e-8;

using namespace Dune;


template<int dim, int order>
bool testShapeFunctionValue(const GeometryType& gt,
                            const FieldVector<double, dim> &pos,
                            int comp, double expected)
{
  MonomialLocalFiniteElement<double,double,dim,order> shapeFunctionSet(gt);
  std::vector<FieldVector<double,1> > out;
  shapeFunctionSet.localBasis().evaluateFunction(pos, out);

  if(std::abs(out[comp][0]-expected) > epsilon) {
    std::cerr << "Bug in shape function of dimension " << dim
              << " and order " << order << " for "
              << gt << "." << std::endl;
    std::cerr << "Value of shape function number " << comp << " at position "
              << pos << " is " << out[comp][0] << " but " << expected
              << " was expected." << std::endl;
    return false;
  }
  return true;
}

int main (int argc, char *argv[])
{
  bool success = true;

  // Do the standard shape function tests
  std::cout << "Monomials are only tested up to order 2 due to the instability of interpolate()." << std::endl;
  Hybrid::forEach(std::make_index_sequence<3>{},[&success](auto i)
  {
    Dune::MonomialLocalFiniteElement<double,double,1,i> monom1d(GeometryTypes::line);
    TEST_FE(monom1d);

    Dune::MonomialLocalFiniteElement<double,double,2,i> monom2d(GeometryTypes::triangle);
    TEST_FE(monom2d);

    Dune::MonomialLocalFiniteElement<double,double,3,i> monom3d(GeometryTypes::tetrahedron);
    TEST_FE(monom3d);
  });

  // Test whether the shape function values are correct
  // dim=1
  success &= testShapeFunctionValue<1,2>(GeometryTypes::line, {0}, 0, 1);
  success &= testShapeFunctionValue<1,2>(GeometryTypes::line, {0}, 1, 0);
  success &= testShapeFunctionValue<1,2>(GeometryTypes::line, {0}, 2, 0);

  success &= testShapeFunctionValue<1,2>(GeometryTypes::line, {0.5}, 0, 1);
  success &= testShapeFunctionValue<1,2>(GeometryTypes::line, {0.5}, 1, .5);
  success &= testShapeFunctionValue<1,2>(GeometryTypes::line, {0.5}, 2, .25);

  success &= testShapeFunctionValue<1,2>(GeometryTypes::line, {1}, 0, 1);
  success &= testShapeFunctionValue<1,2>(GeometryTypes::line, {1}, 1, 1);
  success &= testShapeFunctionValue<1,2>(GeometryTypes::line, {1}, 2, 1);

  // dim=2
  success &= testShapeFunctionValue<2,1>(GeometryTypes::quadrilateral, {0,0}, 0, 1);
  success &= testShapeFunctionValue<2,1>(GeometryTypes::quadrilateral, {0,0}, 1, 0);
  success &= testShapeFunctionValue<2,1>(GeometryTypes::quadrilateral, {0,0}, 2, 0);

  success &= testShapeFunctionValue<2,1>(GeometryTypes::quadrilateral, {0.5,0.5}, 0, 1);
  success &= testShapeFunctionValue<2,1>(GeometryTypes::quadrilateral, {0.5,0.5}, 1, .5);
  success &= testShapeFunctionValue<2,1>(GeometryTypes::quadrilateral, {0.5,0.5}, 2, .5);

  success &= testShapeFunctionValue<2,1>(GeometryTypes::quadrilateral, {1,1}, 0, 1);
  success &= testShapeFunctionValue<2,1>(GeometryTypes::quadrilateral, {1,1}, 1, 1);
  success &= testShapeFunctionValue<2,1>(GeometryTypes::quadrilateral, {1,1}, 2, 1);

  return success ? 0 : 1;
}
