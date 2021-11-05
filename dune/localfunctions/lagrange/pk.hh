// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/* vim: set ai expandtab sw=4 ts=4: */
#ifndef DUNE_PK_LOCALFINITEELEMENT_HH
#define DUNE_PK_LOCALFINITEELEMENT_HH

#include<dune/localfunctions/lagrange/lagrangesimplex.hh>

#warning This header is deprecated

namespace Dune
{
  /**
   * \deprecated This class is obsolete. Please use LagrangeSimplexLocalFiniteElement instead!
   */
  template<class D, class R, int d, int k>
  using PkLocalFiniteElement
    [[deprecated("use LagrangeSimplexLocalFiniteElement instead")]]
    = LagrangeSimplexLocalFiniteElement<D, R, d, k>;
}

#endif
