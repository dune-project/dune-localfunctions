// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include "config.h"

#include <dune/localfunctions/nedelec/nedelec1stkindsimplex.hh>

#include <dune/localfunctions/test/test-localfe.hh>

using namespace Dune;

int main(int argc, char** argv)
{
  bool success = true;

  // First order on a triangle
  Nedelec1stKindSimplexLocalFiniteElement<double,double,2,1> nedelecLFEMTriangle1stOrder;
  TEST_FE3(nedelecLFEMTriangle1stOrder, DisableNone, 2);

  for (unsigned int s = 0; s < 8; s++)
  {
    Nedelec1stKindSimplexLocalFiniteElement<double,double,2,1> nedelecLFEMTriangle1stOrder(s);
    TEST_FE3(nedelecLFEMTriangle1stOrder, DisableNone, 2);
  }

  // First order on a tetrahedron
  Nedelec1stKindSimplexLocalFiniteElement<double,double,3,1> nedelecLFEMTetrahedron1stOrder;
  TEST_FE3(nedelecLFEMTetrahedron1stOrder, DisableNone, 2);

  for (unsigned int s = 0; s < 64; s++)
  {
    Nedelec1stKindSimplexLocalFiniteElement<double,double,3,1> nedelecLFEMTetrahedron1stOrder(s);
    TEST_FE3(nedelecLFEMTetrahedron1stOrder, DisableNone, 2);
  }

  return success ? 0 : 1;
}
