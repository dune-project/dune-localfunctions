// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#include "config.h"

#include <dune/localfunctions/nedelec/nedelec1stkindsimplex.hh>
#include <dune/localfunctions/nedelec/nedelec1stkindcube.hh>

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

  // First order on a square
  Nedelec1stKindCubeLocalFiniteElement<double,double,2,1> nedelecLFEMSquare1stOrder;
  TEST_FE3(nedelecLFEMSquare1stOrder, DisableNone, 2);

  for (unsigned int s = 0; s < 16; s++)
  {
    Nedelec1stKindCubeLocalFiniteElement<double,double,2,1> nedelecLFEMSquare1stOrder(s);
    TEST_FE3(nedelecLFEMSquare1stOrder, DisableNone, 2);
  }

  // First order on a cube
  Nedelec1stKindCubeLocalFiniteElement<double,double,3,1> nedelecLFEMCube1stOrder;
  TEST_FE3(nedelecLFEMCube1stOrder, DisableNone, 2);

  for (unsigned int s = 0; s < 4096; s++)
  {
    Nedelec1stKindCubeLocalFiniteElement<double,double,3,1> nedelecLFEMCube1stOrder(s);
    TEST_FE3(nedelecLFEMCube1stOrder, DisableNone, 2);
  }

  return success ? 0 : 1;
}
