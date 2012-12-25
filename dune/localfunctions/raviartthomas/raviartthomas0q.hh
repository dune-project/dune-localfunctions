// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_RAVIARTTHOMAS0QLOCALFINITEELEMENT_HH
#define DUNE_RAVIARTTHOMAS0QLOCALFINITEELEMENT_HH

#warning This header is deprecated, please use\
  dune/localfunctions/raviartthomas/raviartthomascube.hh instead

#include "raviartthomascube.hh"

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
   *
   * \deprecated This class is deprecated and will be removed after Dune 2.3.
   *             Use RaviartThomasCubeLocalFiniteElement<D, R, dim, 0> instead.
   */
  template<class D, class R, int dim>
  class RT0QLocalFiniteElement;

  /**
   * \brief Specialization on 2d quadrilaterals for lowest order Raviart-Thomas shape functions.
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   *
   * \deprecated This class is deprecated and will be removed after Dune 2.3.
   *             Use RaviartThomasCubeLocalFiniteElement<D, R, 2, 0> instead.
   */
  template<class D, class R>
  class
  DUNE_DEPRECATED_MSG("Use RaviartThomasCubeLocalFiniteElement<D, R, 2, 0> instead")
  RT0QLocalFiniteElement<D, R, 2>
    : public RaviartThomasCubeLocalFiniteElement<D, R, 2, 0>
  {
  public:
    RT0QLocalFiniteElement ()
      : RaviartThomasCubeLocalFiniteElement<D, R, 2, 0>::RaviartThomasCubeLocalFiniteElement()
    {}

    RT0QLocalFiniteElement (int s)
      : RaviartThomasCubeLocalFiniteElement<D, R, 2, 0>::RaviartThomasCubeLocalFiniteElement(s)
    {}
  };

  /**
   * \brief Specialization on 3d quadrilaterals for lowest order Raviart-Thomas shape functions.
   *
   * \tparam D Type to represent the field in the domain.
   * \tparam R Type to represent the field in the range.
   *
   * \deprecated This class is deprecated and will be removed after Dune 2.3.
   *             Use RaviartThomasCubeLocalFiniteElement<D, R, 3, 0> instead.
   */
  template<class D, class R>
  class
  DUNE_DEPRECATED_MSG("Use RaviartThomasCubeLocalFiniteElement<D, R, 3, 0> instead")
  RT0QLocalFiniteElement<D, R, 3>
    : public RaviartThomasCubeLocalFiniteElement<D, R, 3, 0>
  {
  public:
    RT0QLocalFiniteElement ()
      : RaviartThomasCubeLocalFiniteElement<D, R, 3, 0>::RaviartThomasCubeLocalFiniteElement()
    {}

    RT0QLocalFiniteElement (int s)
      : RaviartThomasCubeLocalFiniteElement<D, R, 3, 0>::RaviartThomasCubeLocalFiniteElement(s)
    {}
  };

}

#endif
