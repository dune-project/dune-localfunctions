// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_Q2_LOCALFINITEELEMENT_HH
#define DUNE_Q2_LOCALFINITEELEMENT_HH

#include <dune/common/typetraits.hh>
#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localtoglobaladaptors.hh>

#include <dune/localfunctions/lagrange/lagrangecube.hh>
#include <dune/localfunctions/lagrange/lagrangesimplex.hh>

namespace Dune
{
  //! Factory for global-valued Q23D elements
  /**
   * \tparam Geometry Type of the geometry.  Used to extract the domain field
   *                  type.
   * \tparam RF       Range field type.
   */
  template<class Geometry, class RF>
  class Q2FiniteElementFactory :
    public ScalarLocalToGlobalFiniteElementAdaptorFactory<
      typename std::conditional<Geometry::mydimension == 1,
        LagrangeSimplexLocalFiniteElement<typename Geometry::ctype, RF, 1, 2>,
        LagrangeCubeLocalFiniteElement<typename Geometry::ctype, RF, Geometry::mydimension, 2> >::type,
      Geometry>
  {
    typedef typename std::conditional<Geometry::mydimension == 1,
      LagrangeSimplexLocalFiniteElement<typename Geometry::ctype, RF, 1, 2>,
      LagrangeCubeLocalFiniteElement<typename Geometry::ctype, RF, Geometry::mydimension, 2> >::type LFE;
    typedef ScalarLocalToGlobalFiniteElementAdaptorFactory<LFE, Geometry> Base;

    static const LFE lfe;

  public:
    //! default constructor
    Q2FiniteElementFactory() : Base(lfe) {}
  };

  template<class Geometry, class RF>
  const typename Q2FiniteElementFactory<Geometry, RF>::LFE
  Q2FiniteElementFactory<Geometry, RF>::lfe;
}

#endif
