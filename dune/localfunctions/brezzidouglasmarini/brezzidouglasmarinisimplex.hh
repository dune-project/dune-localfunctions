// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI_BREZZIDOUGLASMARINISIMPLEX_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI_BREZZIDOUGLASMARINISIMPLEX_HH

#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini1simplex2d.hh>
#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini2simplex2d.hh>


namespace Dune
{
  /**
   * \brief Brezzi-Douglas-Marini local finite element for simplices
   *
   * \tparam D Number type to represent domain coordinates
   * \tparam R Number type to represent shape function values
   * \tparam dim Dimension of the reference elements, currently only 2 is supported
   * \tparam order Polynomial order of the element
   */
  template<class D, class R, unsigned int dim, unsigned int order>
  class BrezziDouglasMariniSimplexLocalFiniteElement;

  /**
   * \brief Brezzi-Douglas-Marini local finite elements for simplices with dimension 2 and order 1.
   */
  template<class D, class R>
  class BrezziDouglasMariniSimplexLocalFiniteElement<D, R, 2, 1>
    : public BDM1Simplex2DLocalFiniteElement<D, R>
  {
  public:
    /** \brief Default constructor */
    BrezziDouglasMariniSimplexLocalFiniteElement()
    {}

    /**
     * \brief Constructor with a set of edge orientations
     *
     * \param s Bitfield of size 3 giving the orientations of the three element edges
     */
    BrezziDouglasMariniSimplexLocalFiniteElement(int s)
      : BDM1Simplex2DLocalFiniteElement<D, R>::BDM1Simplex2DLocalFiniteElement(s)
    {}
  };

  /**
   * \brief Brezzi-Douglas-Marini local finite elements for simplices with dimension 2 and order 2.
   */
  template<class D, class R>
  class BrezziDouglasMariniSimplexLocalFiniteElement<D, R, 2, 2>
    : public BDM2Simplex2DLocalFiniteElement<D, R>
  {
  public:
    /** \brief Default constructor */
    BrezziDouglasMariniSimplexLocalFiniteElement()
    {}

    /**
     * \brief Constructor with a set of edge orientations
     *
     * \param s Bitfield of size 3 giving the orientations of the three element edges
     */
    BrezziDouglasMariniSimplexLocalFiniteElement(int s)
      : BDM2Simplex2DLocalFiniteElement<D, R>::BDM2Simplex2DLocalFiniteElement(s)
    {}
  };

} // namespace Dune

#endif // #ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI_BREZZIDOUGLASMARINISIMPLEX_HH
