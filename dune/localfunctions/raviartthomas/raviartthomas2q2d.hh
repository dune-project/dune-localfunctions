// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS2Q2DLOCALFINITEELEMENT_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS2Q2DLOCALFINITEELEMENT_HH

#include "raviartthomas2cube2d.hh"

#warning This header is deprecated, please use\
  dune/localfunctions/raviartthomas/raviartthomas2cube2d.hh instead

namespace Dune
{
  /**
   * \brief Second order Raviart-Thomas shape functions on cubes.
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   *
   * \deprecated This class is deprecated and will be removed after Dune 2.3.
   *             Use RT2Cube2DLocalFiniteElement instead.
   */
  template<class D, class R>
  class
  DUNE_DEPRECATED_MSG("Use RT2Cube2DLocalFiniteElement instead")
  RT2Q2DLocalFiniteElement
    : public RT2Cube2DLocalFiniteElement<D, R>
  {};
}
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS2Q2DLOCALFINITEELEMENT_HH
