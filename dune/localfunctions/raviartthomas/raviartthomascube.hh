// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS_CUBE_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS_CUBE_HH

#include "raviartthomas0cube2d.hh"
#include "raviartthomas0cube3d.hh"
#include "raviartthomas1cube2d.hh"
#include "raviartthomas1cube3d.hh"
#include "raviartthomas2cube2d.hh"
#include "raviartthomas3cube2d.hh"
#include "raviartthomas4cube2d.hh"

/**
 * \file
 * \brief Convenience header that includes all available Raviart-Thomas
 *        local finite elements for cubes.
 */

namespace Dune
{
  /**
   * \brief Raviart-Thomas local finite elements for cubes.
   *
   * Convenience class to access all implemented Raviart-Thomas local
   * finite elements for cubes.
   * Generic Raviart-Thomas local finite elements for cubes like for
   * simpleces could be written.
   *
   * \ingroup RaviartThomas
   *
   * \tparam D type to represent the field in the domain.
   * \tparam R type to represent the field in the range.
   * \tparam dim dimension of the reference elements, must be 2 or 3.
   * \tparam order order of the element, depending on \a dim it can be 0, 1, or 2.
   */
  template<class D, class R, unsigned int dim, unsigned int order>
  class RaviartThomasCubeLocalFiniteElement;

  /**
   * \brief Raviart-Thomas local finite elements for cubes with dimension 2 and order 0.
   */
  template<class D, class R>
  class RaviartThomasCubeLocalFiniteElement<D, R, 2, 0>
    : public RT0Cube2DLocalFiniteElement<D, R>
  {
  public:
    RaviartThomasCubeLocalFiniteElement()
      : RT0Cube2DLocalFiniteElement<D, R>::RT0Cube2DLocalFiniteElement()
    {}

    RaviartThomasCubeLocalFiniteElement(int s)
      : RT0Cube2DLocalFiniteElement<D, R>::RT0Cube2DLocalFiniteElement(s)
    {}
  };

  /**
   * \brief Raviart-Thomas local finite elements for cubes with dimension 2 and order 1.
   */
  template<class D, class R>
  class RaviartThomasCubeLocalFiniteElement<D, R, 2, 1>
    : public RT1Cube2DLocalFiniteElement<D, R>
  {
  public:
    RaviartThomasCubeLocalFiniteElement()
      : RT1Cube2DLocalFiniteElement<D, R>::RT1Cube2DLocalFiniteElement()
    {}

    RaviartThomasCubeLocalFiniteElement(int s)
      : RT1Cube2DLocalFiniteElement<D, R>::RT1Cube2DLocalFiniteElement(s)
    {}
  };

  /**
   * \brief Raviart-Thomas local finite elements for cubes with dimension 2 and order 2.
   */
  template<class D, class R>
  class RaviartThomasCubeLocalFiniteElement<D, R, 2, 2>
    : public RT2Cube2DLocalFiniteElement<D, R>
  {
  public:
    RaviartThomasCubeLocalFiniteElement()
      : RT2Cube2DLocalFiniteElement<D, R>::RT2Cube2DLocalFiniteElement()
    {}

    RaviartThomasCubeLocalFiniteElement(int s)
      : RT2Cube2DLocalFiniteElement<D, R>::RT2Cube2DLocalFiniteElement(s)
    {}
  };

  /**
   * \brief Raviart-Thomas local finite elements for cubes with dimension 2 and order 3.
   */
  template<class D, class R>
  class RaviartThomasCubeLocalFiniteElement<D, R, 2, 3>
    : public RT3Cube2DLocalFiniteElement<D, R>
  {
  public:
    RaviartThomasCubeLocalFiniteElement()
      : RT3Cube2DLocalFiniteElement<D, R>::RT3Cube2DLocalFiniteElement()
    {}

    RaviartThomasCubeLocalFiniteElement(int s)
      : RT3Cube2DLocalFiniteElement<D, R>::RT3Cube2DLocalFiniteElement(s)
    {}
  };

  /**
   * \brief Raviart-Thomas local finite elements for cubes with dimension 2 and order 4.
   */
  template<class D, class R>
  class RaviartThomasCubeLocalFiniteElement<D, R, 2, 4>
    : public RT4Cube2DLocalFiniteElement<D, R>
  {
  public:
    RaviartThomasCubeLocalFiniteElement()
      : RT4Cube2DLocalFiniteElement<D, R>::RT4Cube2DLocalFiniteElement()
    {}

    RaviartThomasCubeLocalFiniteElement(int s)
      : RT4Cube2DLocalFiniteElement<D, R>::RT4Cube2DLocalFiniteElement(s)
    {}
  };

  /**
   * \brief Raviart-Thomas local finite elements for cubes with dimension 3 and order 0.
   */
  template<class D, class R>
  class RaviartThomasCubeLocalFiniteElement<D, R, 3, 0>
    : public RT0Cube3DLocalFiniteElement<D, R>
  {
  public:
    RaviartThomasCubeLocalFiniteElement()
      : RT0Cube3DLocalFiniteElement<D, R>::RT0Cube3DLocalFiniteElement()
    {}

    RaviartThomasCubeLocalFiniteElement(int s)
      : RT0Cube3DLocalFiniteElement<D, R>::RT0Cube3DLocalFiniteElement(s)
    {}
  };

  /**
   * \brief Raviart-Thomas local finite elements for cubes with dimension 3 and order 1.
   */
  template<class D, class R>
  class RaviartThomasCubeLocalFiniteElement<D, R, 3, 1>
    : public RT1Cube3DLocalFiniteElement<D, R>
  {
  public:
    RaviartThomasCubeLocalFiniteElement()
      : RT1Cube3DLocalFiniteElement<D, R>::RT1Cube3DLocalFiniteElement()
    {}

    RaviartThomasCubeLocalFiniteElement(int s)
      : RT1Cube3DLocalFiniteElement<D, R>::RT1Cube3DLocalFiniteElement(s)
    {}
  };
} // namespace Dune

#endif // #ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS_CUBE_HH
