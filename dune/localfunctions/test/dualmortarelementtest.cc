// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#include "config.h"

#include <dune/localfunctions/dualmortarbasis.hh>

#include "test-localfe.hh"

int main(int argc, char** argv)
{
  bool success = true;

  Dune::DualP1LocalFiniteElement<double,double,1> dualp11dlfem;
  TEST_FE(dualp11dlfem);

  Dune::DualP1LocalFiniteElement<double,double,2> dualp12dlfem;
  TEST_FE(dualp12dlfem);

  Dune::DualP1LocalFiniteElement<double,double,3> dualp13dlfem;
  TEST_FE(dualp13dlfem);

  Dune::DualQ1LocalFiniteElement<double,double,1> dualq11dlfem;
  TEST_FE(dualq11dlfem);

  Dune::DualQ1LocalFiniteElement<double,double,2> dualq12dlfem;
  TEST_FE(dualq12dlfem);

  Dune::DualQ1LocalFiniteElement<double,double,3> dualq13dlfem;
  TEST_FE(dualq13dlfem);

  return success ? 0 : 1;
}
