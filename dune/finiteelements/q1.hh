// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_Q1LOCALFINITEELEMENT_HH
#define DUNE_Q1LOCALFINITEELEMENT_HH

#include "p11d.hh"
#include "q12d.hh"
#include "q13d.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R, int d>
  class Q1LocalFiniteElement;

  /** \todo Please doc me !
   */
  template<class D, class R>
  class Q1LocalFiniteElement<D, R, 1>
    : public P11DLocalFiniteElement<D, R>
  {};

  /** \todo Please doc me !
   */
  template<class D, class R>
  class Q1LocalFiniteElement<D, R, 2>
    : public Q12DLocalFiniteElement<D, R>
  {};

  /** \todo Please doc me !
   */
  template<class D, class R>
  class Q1LocalFiniteElement<D, R, 3>
    : public Q13DLocalFiniteElement<D, R>
  {};

}

#endif
