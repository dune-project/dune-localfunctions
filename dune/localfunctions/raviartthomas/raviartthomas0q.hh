// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RAVIARTTHOMAS0QLOCALFINITEELEMENT_HH
#define DUNE_RAVIARTTHOMAS0QLOCALFINITEELEMENT_HH

#include "raviartthomas0q2d.hh"
#include "raviartthomas0q3d.hh"

namespace Dune
{

  /** \todo Please doc me !
   */
  template<class D, class R, int d>
  class RT0QLocalFiniteElement;

  /** \todo Please doc me !
   */
  template<class D, class R>
  class RT0QLocalFiniteElement<D, R, 2>
    : public RT0Q2DLocalFiniteElement<D, R>
  {
  public:
    RT0QLocalFiniteElement () : RT0Q2DLocalFiniteElement<D, R>()
    {}

    RT0QLocalFiniteElement (int s) : RT0Q2DLocalFiniteElement<D, R>(s)
    {}
  };

  /** \todo Please doc me !
   */
  template<class D, class R>
  class RT0QLocalFiniteElement<D, R, 3>
    : public RT0Q3DLocalFiniteElement<D, R>
  {
  public:
    RT0QLocalFiniteElement () : RT0Q3DLocalFiniteElement<D, R>()
    {}

    RT0QLocalFiniteElement (int s) : RT0Q3DLocalFiniteElement<D, R>(s)
    {}
  };

}

#endif
