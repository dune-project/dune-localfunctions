// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RAVIARTTHOMAS0Q3DLOCALFINITEELEMENT_HH
#define DUNE_RAVIARTTHOMAS0Q3DLOCALFINITEELEMENT_HH

#include "raviartthomas0cube3d.hh"

#warning This header is deprecated, please use\
  dune/localfunctions/raviartthomas/raviartthomas0cube3d.hh instead

namespace Dune
{
  /**
   * \deprecated This class is deprecated and will be removed after Dune 2.3.
   *             Use RT0Cube3DLocalFiniteElement instead.
   */
  template<class D, class R>
  class
  DUNE_DEPRECATED_MSG("Use RT0Cube3DLocalFiniteElement instead")
  RT0Q3DLocalFiniteElement
    : public RT0Cube3DLocalFiniteElement<D, R>
  {};
}
#endif
