// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#include "config.h"

#include <iostream>

#include <dune/localfunctions/brezzidouglasfortinmarini/bdfmcube.hh>

#include <dune/localfunctions/test/test-localfe.hh>

int main(int argc, char** argv)
{
  bool success = true;

  Dune::BDFMCubeLocalFiniteElement<double,double, 2, 1> bdfm1cube2dlfem(1);
  TEST_FE(bdfm1cube2dlfem);

  Dune::BDFMCubeLocalFiniteElement<double,double, 2, 2> bdfm2cube2dlfem(1);
  TEST_FE(bdfm2cube2dlfem);

  Dune::BDFMCubeLocalFiniteElement<double,double, 2, 3> bdfm3cube2dlfem(1);
  TEST_FE(bdfm3cube2dlfem);

  return success ? 0 : 1;
}
