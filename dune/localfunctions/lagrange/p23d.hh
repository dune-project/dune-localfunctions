// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P2_3DLOCALFINITEELEMENT_HH
#define DUNE_P2_3DLOCALFINITEELEMENT_HH

#include <dune/localfunctions/lagrange/lagrangesimplex.hh>

namespace Dune
{

  /** \brief Second-order Lagrange local finite element on the reference tetrahedron
   *
   * \tparam D Number type used for domain coordinates
   * \tparam R Number type used for shape function values
   *
   * \deprecated This class is obsolete. Please use LagrangeSimplexLocalFiniteElement instead!
   */
  template<class D, class R>
  using P23DLocalFiniteElement = LagrangeSimplexLocalFiniteElement<D,R,3,2>;

}

#endif
