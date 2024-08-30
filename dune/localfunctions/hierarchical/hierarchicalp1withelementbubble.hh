// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_HIERARCHICAL_P1_WITH_ELEMENTBUBBLE_LOCALFINITEELEMENT_HH
#define DUNE_HIERARCHICAL_P1_WITH_ELEMENTBUBBLE_LOCALFINITEELEMENT_HH

#include <dune/localfunctions/enriched/simplexp1bubble.hh>

namespace Dune
{
  /**
   * \brief Linear Lagrange functions enriched with an element bubble function.
   * \copydoc SimplexP1BubbleLocalFiniteElement
   */
  template<class D, class R, int dim>
  using HierarchicalP1WithElementBubbleLocalFiniteElement
    = SimplexP1BubbleLocalFiniteElement<D,R,dim>;
}

#endif // DUNE_HIERARCHICAL_P1_WITH_ELEMENTBUBBLE_LOCALFINITEELEMENT_HH
