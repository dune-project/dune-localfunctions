// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGELFECACHE_HH
#define DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGELFECACHE_HH

#include <tuple>
#include <utility>

#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

#include <dune/localfunctions/lagrange/lagrangecube.hh>
#include <dune/localfunctions/lagrange/lagrangeprism.hh>
#include <dune/localfunctions/lagrange/lagrangepyramid.hh>
#include <dune/localfunctions/lagrange/lagrangesimplex.hh>
#include <dune/localfunctions/lagrange/p0.hh>
#include <dune/localfunctions/common/localfiniteelementvariantcache.hh>


namespace Dune {



namespace Impl {

  // Provide implemented Lagrange local finite elements

  template<class D, class R, std::size_t dim, std::size_t order>
  struct ImplementedLagrangeFiniteElements : public FixedDimLocalGeometryTypeIndex<dim>
  {
    using FixedDimLocalGeometryTypeIndex<dim>::index;
    static auto getImplementations()
    {
      return std::make_tuple(
        std::make_pair(index(GeometryTypes::simplex(dim)), []() { return LagrangeSimplexLocalFiniteElement<D,R,dim,order>(); }),
        std::make_pair(index(GeometryTypes::cube(dim)),    []() { return LagrangeCubeLocalFiniteElement<D,R,dim,order>(); })
      );
    }
  };

  template<class D, class R, std::size_t dim>
  struct ImplementedLagrangeFiniteElements<D,R,dim,0> : public FixedDimLocalGeometryTypeIndex<dim>
  {
    using FixedDimLocalGeometryTypeIndex<dim>::index;
    static auto getImplementations()
    {
      return std::make_tuple(
        std::make_pair(index(GeometryTypes::simplex(dim)), []() { return P0LocalFiniteElement<D,R,dim>(GeometryTypes::simplex(dim)); }),
        std::make_pair(index(GeometryTypes::cube(dim)),    []() { return P0LocalFiniteElement<D,R,dim>(GeometryTypes::cube(dim)); }),
        std::make_pair(index(GeometryTypes::none(dim)),    []() { return P0LocalFiniteElement<D,R,dim>(GeometryTypes::none(dim)); })
      );
    }
  };

  template<class D, class R>
  struct ImplementedLagrangeFiniteElements<D,R,3,0> : public FixedDimLocalGeometryTypeIndex<3>
  {
    using FixedDimLocalGeometryTypeIndex<3>::index;
    static auto getImplementations()
    {
      return std::make_tuple(
        std::make_pair(index(GeometryTypes::tetrahedron), []() { return P0LocalFiniteElement<D,R,3>(GeometryTypes::tetrahedron); }),
        std::make_pair(index(GeometryTypes::hexahedron),  []() { return P0LocalFiniteElement<D,R,3>(GeometryTypes::hexahedron); }),
        std::make_pair(index(GeometryTypes::prism),       []() { return P0LocalFiniteElement<D,R,3>(GeometryTypes::prism); }),
        std::make_pair(index(GeometryTypes::pyramid),     []() { return P0LocalFiniteElement<D,R,3>(GeometryTypes::pyramid); })
      );
    }
  };

  template<class D, class R>
  struct ImplementedLagrangeFiniteElements<D,R,3,1> : public FixedDimLocalGeometryTypeIndex<3>
  {
    using FixedDimLocalGeometryTypeIndex<3>::index;
    static auto getImplementations()
    {
      return std::make_tuple(
        std::make_pair(index(GeometryTypes::tetrahedron), []() { return LagrangeSimplexLocalFiniteElement<D,R,3,1>(); }),
        std::make_pair(index(GeometryTypes::hexahedron),  []() { return LagrangeCubeLocalFiniteElement<D,R,3,1>(); }),
        std::make_pair(index(GeometryTypes::prism),       []() { return LagrangePrismLocalFiniteElement<D,R,1>(); }),
        std::make_pair(index(GeometryTypes::pyramid),     []() { return LagrangePyramidLocalFiniteElement<D,R,1>(); })
      );
    }
  };

  template<class D, class R>
  struct ImplementedLagrangeFiniteElements<D,R,3,2> : public FixedDimLocalGeometryTypeIndex<3>
  {
    using FixedDimLocalGeometryTypeIndex<3>::index;
    static auto getImplementations()
    {
      return std::make_tuple(
        std::make_pair(index(GeometryTypes::tetrahedron), []() { return LagrangeSimplexLocalFiniteElement<D,R,3,2>(); }),
        std::make_pair(index(GeometryTypes::hexahedron),  []() { return LagrangeCubeLocalFiniteElement<D,R,3,2>(); }),
        std::make_pair(index(GeometryTypes::prism),       []() { return LagrangePrismLocalFiniteElement<D,R,2>(); }),
        std::make_pair(index(GeometryTypes::pyramid),     []() { return LagrangePyramidLocalFiniteElement<D,R,2>(); })
      );
    }
  };

} // namespace Impl



/** \brief A cache that stores all available Pk/Qk like local finite elements for the given dimension and order
 *
 * An interface for dealing with different vertex orders is currently missing.
 *
 * \tparam D Type used for domain coordinates
 * \tparam R Type used for shape function values
 * \tparam dim Element dimension
 * \tparam order Element order
 *
 * The cached finite element implementations can be obtained using get(GeometryType).
 */
template<class D, class R, std::size_t dim, std::size_t order>
using LagrangeLocalFiniteElementCache = LocalFiniteElementVariantCache<Impl::ImplementedLagrangeFiniteElements<D,R,dim,order>>;



} // namespace Dune




#endif // DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGELFECACHE_HH
