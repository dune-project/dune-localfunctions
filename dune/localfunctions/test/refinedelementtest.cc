// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#include "config.h"

#include <iostream>

#include <dune/localfunctions/refined/refinedp1.hh>
#include <dune/localfunctions/refined/refinedp0.hh>

#include <dune/localfunctions/test/test-localfe.hh>

int main(int argc, char** argv)
{
  bool success = true;

  Dune::RefinedP1LocalFiniteElement<double,double,1> refp11dlfem;
  TEST_FE(refp11dlfem);

  Dune::RefinedP1LocalFiniteElement<double,double,2> refp12dlfem;
  TEST_FE(refp12dlfem);

  Dune::RefinedP1LocalFiniteElement<double,double,3> refp13dlfem;
  TEST_FE(refp13dlfem);

  Dune::RefinedP0LocalFiniteElement<double,double,1> refp01dlfem;
  TEST_FE(refp01dlfem);

  Dune::RefinedP0LocalFiniteElement<double,double,2> refp02dlfem;
  TEST_FE(refp02dlfem);

  return success ? 0 : 1;
}
