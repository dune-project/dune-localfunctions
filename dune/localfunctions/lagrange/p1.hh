// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P1LOCALFINITEELEMENT_HH
#define DUNE_P1LOCALFINITEELEMENT_HH

#include <dune/localfunctions/lagrange/lagrangesimplex.hh>

namespace Dune
{

  /** \brief The local p1 finite element on simplices
      \tparam D Domain data type
      \tparam R Range data type
      \tparam dim Dimension of the simplex

      \deprecated This class is obsolete. Please use LagrangeSimplexLocalFiniteElement instead!
   */
  template<class D, class R, int dim>
  using P1LocalFiniteElement = LagrangeSimplexLocalFiniteElement<D,R,dim,1>;

}

#endif
