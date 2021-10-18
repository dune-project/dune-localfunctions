// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_LOCALFUNCTIONS_QK_LOCALFINITEELEMENT_HH
#define DUNE_LOCALFUNCTIONS_QK_LOCALFINITEELEMENT_HH

#include <dune/localfunctions/lagrange/lagrangecube.hh>

#warning This header is deprecated

namespace Dune
{
  /** \brief General Lagrange finite element for cubes with arbitrary dimension and polynomial order
   *   \note The general class QkLocalCoefficients is available for k>0 in dimensions 2 and 3 only
   *
   * \tparam D type used for domain coordinates
   * \tparam R type used for function values
   * \tparam d dimension of the reference element
   * \tparam k polynomial order
   *
   * \deprecated This class is deprecated!  Please use LagrangeCubeLocalFiniteElement instead.
   */
  template<class D, class R, int d, int k>
  using QkLocalFiniteElement
    [[deprecated("use LagrangeCubeLocalFiniteElement instead")]]
    = LagrangeCubeLocalFiniteElement<D,R,d,k>;

}

#endif
