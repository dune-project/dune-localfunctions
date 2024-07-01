// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#include <dune/localfunctions/enriched.hh>

#include <dune/localfunctions/test/test-localfe.hh>

int main(int argc, char** argv)
{
  bool success = true;

  Dune::SimplexP1BubbleLocalFiniteElement<double,double,1> simplexp1b_dim1;
  TEST_FE(simplexp1b_dim1);

  Dune::SimplexP1BubbleLocalFiniteElement<double,double,2> simplexp1b_dim2;
  TEST_FE(simplexp1b_dim2);

  Dune::SimplexP1BubbleLocalFiniteElement<double,double,3> simplexp1b_dim3;
  TEST_FE(simplexp1b_dim3);

  return success ? 0 : 1;
}
