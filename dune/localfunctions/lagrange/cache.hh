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

  //! Construct an empty cache.
  explicit DynamicLagrangeLocalFiniteElementCache (unsigned int order)
    : order_(order)
    , data_()
  {}

  /**
   * \brief Obtain the cached local finite-element.
   *
   * This function might first construct the local finite-element to the polynomial
   * order specified in the constructor of the cache, if it is not yet cached.
   **/
  const FiniteElementType& get (GeometryType type) const
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
class StaticLagrangeLocalFiniteElementCache
{
  struct UnknownToplogy {};

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

  //! Construct the local-finite element for the order specified as template parameter.
  explicit StaticLagrangeLocalFiniteElementCache (std::integral_constant<std::size_t,order> = {})
  {
    lfe_.emplace();
  }

  //! Obtain the cached local finite-element.
  const FiniteElementType& get ([[maybe_unused]] GeometryType type) const
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
 * \tparam Domain Type used for domain coordinates
 * \tparam Range Type used for shape function values
 * \tparam dim Element dimension
 * \tparam order Element order
 *
 * The cached finite element implementations can be obtained using get(GeometryType).
 *
 * \note This is a specialization of the fixed-geometry type LFE cache for the ID `GeometryType::Id(~0u)`.
 * This is given by the default `topologyId` in the capability `Dune::Capabilities::hasSingleGeometryType`
 * that can be extracted from grids with support for mixed geometry types.
 */
template <class Domain, class Range, std::size_t dim, std::size_t order>
class StaticLagrangeLocalFiniteElementCache<GeometryType::Id(~0u), Domain, Range, dim, order>
    : public LagrangeLocalFiniteElementCache<Domain,Range,dim,order>
{
  using Base = LagrangeLocalFiniteElementCache<Domain,Range,dim,order>;
public:
  using Base::Base;
};

} // namespace Dune

#endif // DUNE_LOCALFUNCTIONS_LAGRANGE_CACHE_HH
