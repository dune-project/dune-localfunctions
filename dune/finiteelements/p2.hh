// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P2_LOCALFINITEELEMENT_HH
#define DUNE_P2_LOCALFINITEELEMENT_HH

#include "pk2d.hh"
#include "p23d.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R, int d>
  class P2LocalFiniteElement;

  /** \todo Please doc me !
   */
  template<class D, class R>
  class P2LocalFiniteElement<D, R, 2>
    : public Pk2DLocalFiniteElement<D, R, 2>
  {};

  /** \todo Please doc me !
   */
  template<class D, class R>
  class P2LocalFiniteElement<D, R, 3>
    : public P23DLocalFiniteElement<D, R>
  {};

}

#endif
