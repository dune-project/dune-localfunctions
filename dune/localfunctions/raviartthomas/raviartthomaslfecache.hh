// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS_RAVIARTTHOMASLFECACHE_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS_RAVIARTTHOMASLFECACHE_HH

#include <tuple>
#include <utility>

#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

#include <dune/localfunctions/raviartthomas.hh>
#include <dune/localfunctions/common/localfiniteelementvariantcache.hh>

namespace Dune {

namespace Impl {

  // Provide implemented Raviart-Thomas local finite elements

  template<class D, class R, std::size_t dim, std::size_t order>
  struct ImplementedRaviartThomasLocalFiniteElements
  {};

  template<class D, class R>
  struct ImplementedRaviartThomasLocalFiniteElements<D,R,2,0> : public FixedDimLocalGeometryTypeIndex<2>
  {
    using FixedDimLocalGeometryTypeIndex<2>::index;
    static auto getImplementations()
    {
      return std::make_tuple(
        std::make_pair(index(GeometryTypes::triangle),      []() { return RT02DLocalFiniteElement<D,R>(); }),
        std::make_pair(index(GeometryTypes::quadrilateral), []() { return RT0Cube2DLocalFiniteElement<D,R>(); })
      );
    }
  };

  template<class D, class R>
  struct ImplementedRaviartThomasLocalFiniteElements<D,R,2,1> : public FixedDimLocalGeometryTypeIndex<2>
  {
    using FixedDimLocalGeometryTypeIndex<2>::index;
    static auto getImplementations()
    {
      return std::make_tuple(
        std::make_pair(index(GeometryTypes::triangle),      []() { return RT12DLocalFiniteElement<D,R>(); }),
        std::make_pair(index(GeometryTypes::quadrilateral), []() { return RT1Cube2DLocalFiniteElement<D,R>(); })
      );
    }
  };

  template<class D, class R>
  struct ImplementedRaviartThomasLocalFiniteElements<D,R,2,2> : public FixedDimLocalGeometryTypeIndex<2>
  {
    using FixedDimLocalGeometryTypeIndex<2>::index;
    static auto getImplementations()
    {
      return std::make_tuple(
        std::make_pair(index(GeometryTypes::quadrilateral), []() { return RT2Cube2DLocalFiniteElement<D,R>(); })
      );
    }
  };

  template<class D, class R>
  struct ImplementedRaviartThomasLocalFiniteElements<D,R,3,0> : public FixedDimLocalGeometryTypeIndex<3>
  {
    using FixedDimLocalGeometryTypeIndex<3>::index;
    static auto getImplementations()
    {
      return std::make_tuple(
        std::make_pair(index(GeometryTypes::tetrahedron), []() { return RT03DLocalFiniteElement<D,R>(); }),
        std::make_pair(index(GeometryTypes::hexahedron), []() { return RT0Cube3DLocalFiniteElement<D,R>(); })
      );
    }
  };

  template<class D, class R>
  struct ImplementedRaviartThomasLocalFiniteElements<D,R,3,1> : public FixedDimLocalGeometryTypeIndex<3>
  {
    using FixedDimLocalGeometryTypeIndex<3>::index;
    static auto getImplementations()
    {
      return std::make_tuple(
        std::make_pair(index(GeometryTypes::hexahedron), []() { RT1Cube3DLocalFiniteElement<D,R>(); })
      );
    }
  };

} // namespace Impl



/** \brief A cache that stores all available Raviart-Thomas local finite elements for the given dimension and order
 *
 * \tparam D Type used for domain coordinates
 * \tparam R Type used for shape function values
 * \tparam dim Element dimension
 * \tparam order Element order
 *
 * The cached finite element implementations can be obtained using get(GeometryType).
 */
template<class D, class R, std::size_t dim, std::size_t order>
using RaviartThomasLocalFiniteElementCache = LocalFiniteElementVariantCache<Impl::ImplementedRaviartThomasLocalFiniteElements<D,R,dim,order>>;

} // namespace Dune

#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS_RAVIARTTHOMASLFECACHE_HH
