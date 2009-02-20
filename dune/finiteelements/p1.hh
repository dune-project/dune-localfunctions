// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P1LOCALFINITEELEMENT_HH
#define DUNE_P1LOCALFINITEELEMENT_HH

#include "p12d.hh"
#include "p13d.hh"

namespace Dune
{

  template<class D, class R, int d>
  class P1LocalFiniteElement;

  template<class D, class R>
  class P1LocalFiniteElement<D, R, 2>
    : public P12DLocalFiniteElement<D, R>
  {};

  template<class D, class R>
  class P1LocalFiniteElement<D, R, 3>
    : public P13DLocalFiniteElement<D, R>
  {};

}

#endif
