// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RAVIARTTHOMAS0QLOCALFINITEELEMENT_HH
#define DUNE_RAVIARTTHOMAS0QLOCALFINITEELEMENT_HH

#include "raviartthomas0cube2d.hh"
#include "raviartthomas0cube3d.hh"

namespace Dune
{

  /**
   * \brief Lowest order Raviart-Thomas shape functions on quadrilaterals.
   *
   * The dimensions d=2 and dim=3 are supported. This is a convenience class
   * to include the Raviart-Thomas-0 basis elements on quadrilaterals.
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   * \tparam dim Dimension.
   */
  template<class D, class R, int dim>
  class RT0CubeLocalFiniteElement;

  /**
   * \brief Specialization on 2d quadrilaterals for lowest order Raviart-Thomas shape functions.
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   */
  template<class D, class R>
  class RT0CubeLocalFiniteElement<D, R, 2>
    : public RT0Cube2DLocalFiniteElement<D, R>
  {
  public:
    RT0CubeLocalFiniteElement () : RT0Cube2DLocalFiniteElement<D, R>()
    {}

    RT0CubeLocalFiniteElement (int s) : RT0Cube2DLocalFiniteElement<D, R>(s)
    {}
  };

  /**
   * \brief Specialization on 3d quadrilaterals for lowest order Raviart-Thomas shape functions.
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   */
  template<class D, class R>
  class RT0CubeLocalFiniteElement<D, R, 3>
    : public RT0Cube3DLocalFiniteElement<D, R>
  {
  public:
    RT0CubeLocalFiniteElement () : RT0Cube3DLocalFiniteElement<D, R>()
    {}

    RT0CubeLocalFiniteElement (int s) : RT0Cube3DLocalFiniteElement<D, R>(s)
    {}
  };

}

#endif
