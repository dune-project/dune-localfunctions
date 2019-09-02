// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PYRAMID_P1_LOCALFINITEELEMENT_HH
#define DUNE_PYRAMID_P1_LOCALFINITEELEMENT_HH

#include <dune/localfunctions/lagrange/lagrangepyramid.hh>

namespace Dune
{

  /** \brief First-order Lagrangian finite element on a three-dimensional pyramid
   *
   * \deprecated Please use LagrangePyramidLocalFiniteElement<D,R,1> instead!
   */
  template<class D, class R>
  using PyramidP1LocalFiniteElement = LagrangePyramidLocalFiniteElement<D,R,1>;

}

#endif
