// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PRISM2_3DLOCALFINITEELEMENT_HH
#define DUNE_PRISM2_3DLOCALFINITEELEMENT_HH

#include <dune/localfunctions/lagrange/lagrangeprism.hh>

#warning This header is deprecated

namespace Dune
{

  /** \brief Second-order Lagrange finite element on a three-dimensional prism
   *
   * \deprecated Please use LagrangePrismLocalFiniteElement<D,R,2> instead!
   */
  template<class D, class R>
  using PrismP2LocalFiniteElement
    [[deprecated("use LagrangePrismLocalFiniteElement instead")]]
    = LagrangePrismLocalFiniteElement<D,R,2>;

}

#endif
