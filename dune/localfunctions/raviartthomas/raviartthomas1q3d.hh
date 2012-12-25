// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS1Q3DLOCALFINITEELEMENT_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS1Q3DLOCALFINITEELEMENT_HH

#include "raviartthomas1cube3d.hh"

#warning This header is deprecated, please use\
  dune/localfunctions/raviartthomas/raviartthomas0cube3d.hh instead

namespace Dune
{
  /**
   * \brief First order Raviart-Thomas shape functions on cubes.
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   *
   * \deprecated This class is deprecated and will be removed after Dune 2.3.
   *             Use RT1Cube3DLocalFiniteElement instead.
   */
  template<class D, class R>
  class
  DUNE_DEPRECATED_MSG("Use RT1Cube3DLocalFiniteElement instead")
  RT1Q3DLocalFiniteElement
    : public RT1Cube3DLocalFiniteElement<D, R>
  {};
}
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS1Q3DLOCALFINITEELEMENT_HH
