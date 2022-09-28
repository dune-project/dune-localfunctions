// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_Pk3DLOCALFINITEELEMENT_HH
#define DUNE_Pk3DLOCALFINITEELEMENT_HH

#include <dune/localfunctions/lagrange/lagrangesimplex.hh>

#warning This header is deprecated

namespace Dune
{

  /** \todo Please doc me !

      \deprecated This class is obsolete. Please use LagrangeSimplexLocalFiniteElement instead!
   */
  template<class D, class R, unsigned int k>
  using Pk3DLocalFiniteElement
    [[deprecated("use LagrangeSimplexLocalFiniteElement instead")]]
    = LagrangeSimplexLocalFiniteElement<D,R,3,k>;

}

#endif
