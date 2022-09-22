// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_PYRAMIDP2_3DLOCALFINITEELEMENT_HH
#define DUNE_PYRAMIDP2_3DLOCALFINITEELEMENT_HH

#include <dune/localfunctions/lagrange/lagrangepyramid.hh>

#warning This header is deprecated

namespace Dune
{

  /** \brief Second-order Lagrangian finite element on a three-dimensional pyramid
   *
   * \deprecated Please use LagrangePyramidLocalFiniteElement<D,R,2> instead!
   */
  template<class D, class R>
  using PyramidP2LocalFiniteElement
    [[deprecated("use LagrangePyramidLocalFiniteElement instead")]]
    = LagrangePyramidLocalFiniteElement<D,R,2>;

}

#endif
