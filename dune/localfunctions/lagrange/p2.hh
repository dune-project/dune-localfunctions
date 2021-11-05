// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P2_LOCALFINITEELEMENT_HH
#define DUNE_P2_LOCALFINITEELEMENT_HH

#include <dune/localfunctions/lagrange/lagrangesimplex.hh>

#warning This header is deprecated

namespace Dune
{

  /** \brief Second-order Lagrange finite element on the reference simplex with compile-time dimension

      \deprecated This class is obsolete. Please use LagrangeSimplexLocalFiniteElement instead!
   */
  template<class D, class R, int d>
  using P2LocalFiniteElement
    [[deprecated("use LagrangeSimplexLocalFiniteElement instead")]]
    = LagrangeSimplexLocalFiniteElement<D,R,d,2>;

}

#endif
