// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#include "config.h"

#include <iostream>

#include <dune/localfunctions/brezzidouglasmarini.hh>

#include <dune/localfunctions/test/test-localfe.hh>

int main(int argc, char** argv)
{
  bool success = true;

  Dune::BrezziDouglasMariniCubeLocalFiniteElement<double,double,2,1> bdm1cube2dlfem(1);
  TEST_FE(bdm1cube2dlfem);

  Dune::BrezziDouglasMariniCubeLocalFiniteElement<double,double,3,1> bdm1cube3dlfem(1);
  // \todo Implement the missing LocalInterpolation
  // DisableRepresentConstants is only set because the test also uses DisableLocalInterpolation internally.
  TEST_FE2(bdm1cube3dlfem, DisableLocalInterpolation + DisableRepresentConstants);

  Dune::BrezziDouglasMariniCubeLocalFiniteElement<double,double,2,2> bdm2cube2dlfem(1);
  TEST_FE(bdm2cube2dlfem);

  Dune::BrezziDouglasMariniSimplexLocalFiniteElement<double,double,2,1> bdm1simplex2dlfem(1);
  TEST_FE(bdm1simplex2dlfem);

  Dune::BrezziDouglasMariniSimplexLocalFiniteElement<double,double,2,2> bdm2simplex2dlfem(1);
  TEST_FE(bdm2simplex2dlfem);

  return success ? 0 : 1;
}
