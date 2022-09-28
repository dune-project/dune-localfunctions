// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI_BREZZIDOUGLASMARINICUBE_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI_BREZZIDOUGLASMARINICUBE_HH

#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini1cube2d.hh>
#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini1cube3d.hh>
#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini2cube2d.hh>


namespace Dune
{
  /**
   * \brief Brezzi-Douglas-Marini local finite element for cubes
   *
   * \tparam D Number type to represent domain coordinates
   * \tparam R Number type to represent shape function values
   * \tparam dim Dimension of the reference elements, must be 2 or 3
   * \tparam order Polynomial order of the element
   */
  template<class D, class R, unsigned int dim, unsigned int order>
  class BrezziDouglasMariniCubeLocalFiniteElement;

  /**
   * \brief Brezzi-Douglas-Marini local finite elements for cubes with dimension 2 and order 1.
   */
  template<class D, class R>
  class BrezziDouglasMariniCubeLocalFiniteElement<D, R, 2, 1>
    : public BDM1Cube2DLocalFiniteElement<D, R>
  {
  public:
    /** \brief Default constructor */
    BrezziDouglasMariniCubeLocalFiniteElement()
    {}

    /**
     * \brief Constructor with a set of edge orientations
     *
     * \param s Bitfield of size 4 giving the orientations of the four element edges
     */
    BrezziDouglasMariniCubeLocalFiniteElement(int s)
      : BDM1Cube2DLocalFiniteElement<D, R>::BDM1Cube2DLocalFiniteElement(s)
    {}
  };

  /**
   * \brief Brezzi-Douglas-Marini local finite elements for cubes with dimension 2 and order 2.
   */
  template<class D, class R>
  class BrezziDouglasMariniCubeLocalFiniteElement<D, R, 2, 2>
    : public BDM2Cube2DLocalFiniteElement<D, R>
  {
  public:
    /** \brief Default constructor */
    BrezziDouglasMariniCubeLocalFiniteElement()
    {}

    /**
     * \brief Constructor with a set of edge orientations
     *
     * \param s Bitfield of size 4 giving the orientations of the four element edges
     */
    BrezziDouglasMariniCubeLocalFiniteElement(int s)
      : BDM2Cube2DLocalFiniteElement<D, R>::BDM2Cube2DLocalFiniteElement(s)
    {}
  };

  /**
   * \brief Brezzi-Douglas-Marini local finite elements for cubes with dimension 3 and order 1.
   */
  template<class D, class R>
  class BrezziDouglasMariniCubeLocalFiniteElement<D, R, 3, 1>
    : public BDM1Cube3DLocalFiniteElement<D, R>
  {
  public:
    /** \brief Default constructor */
    BrezziDouglasMariniCubeLocalFiniteElement()
    {}

    /**
     * \brief Constructor with a set of edge orientations
     *
     * \param s Bitfield of size 6 giving the orientations of the six element facets
     */
    BrezziDouglasMariniCubeLocalFiniteElement(int s)
      : BDM1Cube3DLocalFiniteElement<D, R>::BDM1Cube3DLocalFiniteElement(s)
    {}
  };

} // namespace Dune

#endif // #ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI_BREZZIDOUGLASMARINICUBE_HH
