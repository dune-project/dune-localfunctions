// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RAVIARTTHOMAS0Q2DLOCALFINITEELEMENT_HH
#define DUNE_RAVIARTTHOMAS0Q2DLOCALFINITEELEMENT_HH

#include "raviartthomas0cube2d.hh"

#warning This header is deprecated, please use\
  dune/localfunctions/raviartthomas/raviartthomas0cube2d.hh instead

namespace Dune
{
  /**
   * \deprecated This class is deprecated and will be removed after Dune 2.3.
   *             Use RT0Cube2DLocalFiniteElement instead.
   */
  template<class D, class R>
  class
  DUNE_DEPRECATED_MSG("Use RT0Cube2DLocalFiniteElement instead")
  RT0Q2DLocalFiniteElement
    : public RT0Cube2DLocalFiniteElement<D, R>
  {};
}
#endif
