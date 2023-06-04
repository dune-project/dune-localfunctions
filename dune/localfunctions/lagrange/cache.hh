// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_LAGRANGE_CACHE_HH
#define DUNE_LOCALFUNCTIONS_LAGRANGE_CACHE_HH

#include <map>
#include <optional>
#include <type_traits>

#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

#include <dune/localfunctions/lagrange.hh>
#include <dune/localfunctions/lagrange/equidistantpoints.hh>
#include <dune/localfunctions/lagrange/lagrangecube.hh>
#include <dune/localfunctions/lagrange/lagrangelfecache.hh>
#include <dune/localfunctions/lagrange/lagrangeprism.hh>
#include <dune/localfunctions/lagrange/lagrangepyramid.hh>
#include <dune/localfunctions/lagrange/lagrangesimplex.hh>


namespace Dune {

/** \brief A cache that stores Lagrange finite elements for the given dimension and order.
 *
 * The cache is based on a runtime-order implementation of Lagrange shape functions subject to a
 * given equidistance points-set. The order is given to the class in the constructor.
 *
 * An interface for dealing with different vertex orders is currently missing.
 *
 * \tparam Domain Type used for domain coordinates
 * \tparam Range Type used for shape function values
 * \tparam dim Element dimension
 *
 * The cached finite element implementations can be obtained using get(GeometryType).
 */
template <class Domain, class Range, int dim>
class DynamicLagrangeLocalFiniteElementCache
{
public:
  using FiniteElementType = LagrangeLocalFiniteElement<EquidistantPointSet, dim, Domain, Range>;

  explicit DynamicLagrangeLocalFiniteElementCache (unsigned int order)
    : order_(order)
    , data_()
  {}

  FiniteElementType const& get (GeometryType type) const
  {
    auto [it,_] = data_.try_emplace(type,type,order_);
    return it->second;
  }

private:
  unsigned int order_;
  mutable std::map<GeometryType, FiniteElementType> data_;
};


/** \brief A cache that stores all available Pk/Qk like local finite elements for the given dimension and order
 * for the case that the GeometryType is fixed and has the given Id.
 *
 * \tparam id  The Id of the fixed GeometryType
 * \tparam Domain Type used for domain coordinates
 * \tparam Range Type used for shape function values
 * \tparam dim Element dimension
 * \tparam order Element order
 *
 * The cached finite element implementations can be obtained using get(GeometryType).
 */
template <GeometryType::Id id, class Domain, class Range, std::size_t dim, std::size_t order>
class FixedGeometryTypeLagrangeLocalFiniteElementCache
{
  struct UnknownToplogy {};
  struct UnsupportedVertexMap {};

  static constexpr bool isSimplex = GeometryType(id).isSimplex();
  static constexpr bool isCube = GeometryType(id).isCube();
  static constexpr bool isPrism = GeometryType(id).isPrism();
  static constexpr bool isPyramid = GeometryType(id).isPyramid();

public:
  using FiniteElementType
    = std::conditional_t<isSimplex, LagrangeSimplexLocalFiniteElement<Domain,Range,dim,order>,
      std::conditional_t<isCube,    LagrangeCubeLocalFiniteElement<Domain,Range,dim,order>,
      std::conditional_t<isPrism,   LagrangePrismLocalFiniteElement<Domain,Range,order>,
      std::conditional_t<isPyramid, LagrangePyramidLocalFiniteElement<Domain,Range,order>, UnknownToplogy> > > >;

  using VertexMap
    = std::conditional_t<isSimplex, std::array<unsigned int, dim+1>, UnsupportedVertexMap>;

  FixedGeometryTypeLagrangeLocalFiniteElementCache ()
  {
    lfe_.emplace();
  }

  template <bool isSimplex_ = isSimplex, std::enable_if_t<isSimplex_, int> = 0>
  explicit FixedGeometryTypeLagrangeLocalFiniteElementCache (const VertexMap& vertexMap)
  {
    lfe_.emplace(vertexMap);
  }

  template <bool isSimplex_ = isSimplex, std::enable_if_t<isSimplex_, int> = 0>
  void updateVertexMap (const VertexMap& vertexMap)
  {
    lfe_.emplace(vertexMap);
  }

  FiniteElementType const& get (GeometryType type) const
  {
    assert(GeometryType::Id(type) == id);
    assert(!!lfe_);
    return *lfe_;
  }

private:
  std::optional<FiniteElementType> lfe_{};
};


/** \brief A cache that stores all available Pk/Qk like local finite elements for the given dimension and order
 *
 * An interface for dealing with different vertex orders is currently missing.
 *
 * \tparam Domain Type used for domain coordinates
 * \tparam Range Type used for shape function values
 * \tparam dim Element dimension
 * \tparam order Element order
 *
 * The cached finite element implementations can be obtained using get(GeometryType).
 */
template<class Domain, class Range, std::size_t dim, std::size_t order>
using StaticLagrangeLocalFiniteElementCache = LagrangeLocalFiniteElementCache<Domain,Range,dim,order>;


} // namespace Dune




#endif // DUNE_LOCALFUNCTIONS_LAGRANGE_CACHE_HH
