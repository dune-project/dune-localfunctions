// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#ifndef DUNE_Q1_LOCALFINITEELEMENT_HH
#define DUNE_Q1_LOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localtoglobaladaptors.hh>
#include <dune/localfunctions/lagrange/lagrangecube.hh>

#warning This header is deprecated

namespace Dune
{

  /** \brief The local Q1 finite element on cubes
      \tparam D Domain data type
      \tparam R Range data type
      \tparam dim Dimension of the cube

      \deprecated This class is deprecated!  Please use LagrangeCubeLocalFiniteElement instead.
   */
  template<class D, class R, int dim>
  using Q1LocalFiniteElement
    [[deprecated("use LagrangeCubeLocalFiniteElement instead")]]
    = LagrangeCubeLocalFiniteElement<D,R,dim,1>;


  //! Factory for global-valued Q1 elements
  /**
   * \tparam Geometry Type of the geometry.  Used to extract the domain field
   *                  type and the dimension.
   * \tparam RF       Range field type.
   */
  template<class Geometry, class RF>
  class Q1FiniteElementFactory :
    public ScalarLocalToGlobalFiniteElementAdaptorFactory<
        LagrangeCubeLocalFiniteElement<
            typename Geometry::ctype, RF, Geometry::mydimension, 1
            >,
        Geometry
        >
  {
    typedef LagrangeCubeLocalFiniteElement<
        typename Geometry::ctype, RF, Geometry::mydimension, 1
        > LFE;
    typedef ScalarLocalToGlobalFiniteElementAdaptorFactory<LFE, Geometry> Base;

    static const LFE lfe;

  public:
    //! default constructor
    Q1FiniteElementFactory() : Base(lfe) {}
  };

  template<class Geometry, class RF>
  const typename Q1FiniteElementFactory<Geometry, RF>::LFE
  Q1FiniteElementFactory<Geometry, RF>::lfe;
}

#endif
