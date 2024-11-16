// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#include <functional>

#include <dune/localfunctions/lagrange/lagrangesimplex.hh>
#include <dune/localfunctions/meta/discontinuous.hh>
#include <dune/localfunctions/test/test-localfe.hh>

int main(int argc, char** argv)
{
  bool success = true;

  auto simplexLFE_dim1 = Dune::LagrangeSimplexLocalFiniteElement<double,double,1,3>{};
  auto simplexLFE_dim2 = Dune::LagrangeSimplexLocalFiniteElement<double,double,2,3>{};
  auto simplexLFE_dim3 = Dune::LagrangeSimplexLocalFiniteElement<double,double,3,3>{};

  // transform the Lagrange LFE into a discontinuous version by associating all basis functions
  // with the element interior.

  auto discSimplexLFE_dim1 = Dune::DiscontinuousLocalFiniteElement{simplexLFE_dim1};
  TEST_FE(discSimplexLFE_dim1);

  auto discSimplexLFE_dim2 = Dune::DiscontinuousLocalFiniteElement{simplexLFE_dim2};
  TEST_FE(discSimplexLFE_dim2);

  auto discSimplexLFE_dim3 = Dune::DiscontinuousLocalFiniteElement{simplexLFE_dim3};
  TEST_FE(discSimplexLFE_dim3);

  // Test passing the LFE as reference_wrapper

  auto discSimplexLFE_dim3_ref = Dune::DiscontinuousLocalFiniteElement{std::ref(simplexLFE_dim3)};
  TEST_FE(discSimplexLFE_dim3_ref);

  auto discSimplexLFE_dim3_cref = Dune::DiscontinuousLocalFiniteElement{std::cref(simplexLFE_dim3)};
  TEST_FE(discSimplexLFE_dim3_cref);

  return success ? 0 : 1;
}
