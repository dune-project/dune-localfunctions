// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_NEDELEC_NEDELECSIMPLEX_NEDELECSIMPLEXBASIS_HH

#define DUNE_LOCALFUNCTIONS_NEDELEC_NEDELECSIMPLEX_NEDELECSIMPLEXBASIS_HH

#include <fstream>
#include <dune/common/exceptions.hh>

#include <dune/localfunctions/utility/defaultbasisfactory.hh>
#include "nedelecsimplexinterpolation.hh"
#include "nedelecsimplexprebasis.hh"

namespace Dune
{
  /*
   * `NedelecPreBasisFactory` provides a basis for the Nedelec function space.
   * `NedelecL2InterpolationFactory` provides the linear functionals.
   *
   * `Defaultbasisfactory::create` first builds the function space and the linear functionals.
   * Then the constructor of `BasisMatrix` gets called. There the matrix
   *
   * \begin{equation}
   *   A_{i,j} := N_j(\phi_i)
   * \end{equation}
   *
   * with linear functionals $N_j$ and basisfunctions $\phi_i$ gets assembled.
   * Then the matrix gets inverted and is then used as a coefficent matrix for the standard monomial basis.
   *
   * For more details on the theory see the first chapter "Construction of Local Finite Element Spaces Using the Generic Reference Elements"
   * of the book "Advances in Dune" by Dedner, Flemisch and Klöfkorn published in 2012.
   */

  template< unsigned int dim, class SF, class CF >
  struct NedelecBasisFactory
    : public DefaultBasisFactory< NedelecPreBasisFactory<dim,CF>,
          NedelecL2InterpolationFactory<dim,CF>,
          dim,dim,SF,CF >
  {};
}

#endif // #ifndef DUNE_LOCALFUNCTIONS_NEDELEC_NEDELECSIMPLEX_NEDELECSIMPLEXBASIS_HH
