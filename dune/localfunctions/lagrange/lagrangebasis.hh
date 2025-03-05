// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LAGRANGEBASIS_HH
#define DUNE_LAGRANGEBASIS_HH

#include <fstream>
#include <dune/common/exceptions.hh>

#include <dune/localfunctions/utility/defaultbasisfactory.hh>
#include <dune/localfunctions/utility/monomialbasis.hh>

#include <dune/localfunctions/lagrange/interpolation.hh>

namespace Dune
{
  //! Factory for Lagrange local basis based on a Lagrange point-set
  /**
   * \tparam LP   Template class defining the points for the lagrange interpolation
   * \tparam dim  Dimension of reference elements
   * \tparam D    Domain field-type of the basis functions
   * \tparam R    Range field-type of the basis functions
   * \tpapam SF   Storage field-type for basis matrix
   * \tparam CF   Compute field-type for basis matrix
   **/
  template< template <class,unsigned int> class LP,
      unsigned int dim, class D, class R,
      class SF=R, class CF=SF >
  struct LagrangeBasisFactory
    : public DefaultBasisFactory< MonomialBasisFactory<dim,CF>,
          LagrangeInterpolationFactory<LP,dim,CF>,
          dim,1,D,R,SF,CF >
  {};

}

#endif // #ifndef DUNE_LAGRANGEBASIS_HH
