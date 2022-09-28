// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#include "config.h"

#include <dune/localfunctions/crouzeixraviart.hh>

#include <dune/localfunctions/test/test-localfe.hh>

int main(int argc, char** argv)
{
  bool success = true;

  Dune::CrouzeixRaviartLocalFiniteElement<double,double,1> crozeixRaviart1dLFE;
  TEST_FE3(crozeixRaviart1dLFE, DisableNone, 2 /* difforder */);

  Dune::CrouzeixRaviartLocalFiniteElement<double,double,2> crozeixRaviart2dLFE;
  TEST_FE3(crozeixRaviart2dLFE, DisableNone, 2 /* difforder */);

  Dune::CrouzeixRaviartLocalFiniteElement<double,double,3> crozeixRaviart3dLFE;
  TEST_FE3(crozeixRaviart3dLFE, DisableNone, 2 /* difforder */);

  Dune::CrouzeixRaviartLocalFiniteElement<double,double,4> crozeixRaviart4dLFE;
  TEST_FE3(crozeixRaviart4dLFE, DisableNone, 2 /* difforder */);

  return success ? 0 : 1;
}
