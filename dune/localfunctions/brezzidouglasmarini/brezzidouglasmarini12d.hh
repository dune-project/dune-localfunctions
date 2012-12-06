// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI12DLOCALFINITEELEMENT_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI12DLOCALFINITEELEMENT_HH

#include "brezzidouglasmarini1simplex2d.hh"

#warning This header is deprecated, please use\
  dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini1simplex2d.hh instead

namespace Dune
{

  /**
   * \brief First order Brezzi-Douglas-Marini shape functions on triangles.
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   *
   * \deprecated This class is deprecated and will be removed after Dune 2.3.
   *             Use BDM1Simplex2DLocalFiniteElement instead.
   */
  template<class D, class R>
  class
  DUNE_DEPRECATED_MSG("Use BDM1Simplex2DLocalFiniteElement instead")
  BDM12DLocalFiniteElement
    : public BDM1Simplex2DLocalFiniteElement<D, R>
  {};
}
#endif // DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI12DLOCALFINITEELEMENT_HH
