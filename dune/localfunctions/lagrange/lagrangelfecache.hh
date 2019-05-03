// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGELFECACHE_HH
#define DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGELFECACHE_HH

#include <vector>
#include <tuple>
#include <utility>
#include <type_traits>

#include <dune/common/std/optional.hh>
#include <dune/common/std/type_traits.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/typelist.hh>
#include <dune/common/hybridutilities.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

#include <dune/localfunctions/lagrange/p0.hh>
#include <dune/localfunctions/lagrange/pk.hh>
#include <dune/localfunctions/lagrange/qk.hh>
#include <dune/localfunctions/lagrange/prismp1.hh>
#include <dune/localfunctions/lagrange/prismp2.hh>
#include <dune/localfunctions/lagrange/pyramidp1.hh>
#include <dune/localfunctions/lagrange/pyramidp2.hh>
#include <dune/localfunctions/common/localfiniteelementvariant.hh>


namespace Dune {



namespace Impl {

  // Provide implemented Lagrange local finite elements

  template<class D, class R, std::size_t dim, std::size_t order>
  struct ImplementedLagrangeFiniteElements
  {
    static auto get()
    {
      return std::make_tuple(
        std::make_pair(GeometryTypes::simplex(dim), []() { return PkLocalFiniteElement<D,R,dim,order>(); }),
        std::make_pair(GeometryTypes::cube(dim), []() { return QkLocalFiniteElement<D,R,dim,order>(); })
      );
    }
  };

  template<class D, class R, std::size_t dim>
  struct ImplementedLagrangeFiniteElements<D,R,dim,0>
  {
    static auto get()
    {
      return std::make_tuple(
        std::make_pair(GeometryTypes::simplex(dim), []() { return P0LocalFiniteElement<D,R,dim>(GeometryTypes::simplex(dim)); }),
        std::make_pair(GeometryTypes::cube(dim), []() { return P0LocalFiniteElement<D,R,dim>(GeometryTypes::cube(dim)); }),
        std::make_pair(GeometryTypes::none(dim), []() { return P0LocalFiniteElement<D,R,dim>(GeometryTypes::none(dim)); })
      );
    }
  };

  template<class D, class R>
  struct ImplementedLagrangeFiniteElements<D,R,3,0>
  {
    static auto get()
    {
      return std::make_tuple(
        std::make_pair(GeometryTypes::tetrahedron, []() { return P0LocalFiniteElement<D,R,3>(GeometryTypes::tetrahedron); }),
        std::make_pair(GeometryTypes::hexahedron, []() { return P0LocalFiniteElement<D,R,3>(GeometryTypes::hexahedron); }),
        std::make_pair(GeometryTypes::prism, []() { return P0LocalFiniteElement<D,R,3>(GeometryTypes::prism); }),
        std::make_pair(GeometryTypes::pyramid, []() { return P0LocalFiniteElement<D,R,3>(GeometryTypes::pyramid); })
      );
    }
  };

  template<class D, class R>
  struct ImplementedLagrangeFiniteElements<D,R,3,1>
  {
    static auto get()
    {
      return std::make_tuple(
        std::make_pair(GeometryTypes::tetrahedron, []() { return PkLocalFiniteElement<D,R,3,1>(); }),
        std::make_pair(GeometryTypes::hexahedron, []() { return QkLocalFiniteElement<D,R,3,1>(); }),
        std::make_pair(GeometryTypes::prism, []() { return PrismP1LocalFiniteElement<D,R>(); }),
        std::make_pair(GeometryTypes::pyramid, []() { return PyramidP1LocalFiniteElement<D,R>(); })
      );
    }
  };

  template<class D, class R>
  struct ImplementedLagrangeFiniteElements<D,R,3,2>
  {
    static auto get()
    {
      return std::make_tuple(
        std::make_pair(GeometryTypes::tetrahedron, []() { return PkLocalFiniteElement<D,R,3,2>(); }),
        std::make_pair(GeometryTypes::hexahedron, []() { return QkLocalFiniteElement<D,R,3,2>(); }),
        std::make_pair(GeometryTypes::prism, []() { return PrismP2LocalFiniteElement<D,R>(); }),
        std::make_pair(GeometryTypes::pyramid, []() { return PyramidP2LocalFiniteElement<D,R>(); })
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
 */
template<class D, class R, std::size_t dim, std::size_t order>
class LagrangeLocalFiniteElementCache
{

  template<class LFEImplTuple>
  struct GenerateLFEVariant;

  template<class... LFEImpl>
  struct GenerateLFEVariant<std::tuple<std::pair<GeometryType, LFEImpl>...>>
  {
    using type = UniqueTypes_t<LocalFiniteElementVariant, decltype(std::declval<LFEImpl>()())...>;
  };

  static auto implementedFiniteElements()
  {
    return Impl::ImplementedLagrangeFiniteElements<D,R,dim, order>::get();
  }

public:

  using FiniteElementType = typename GenerateLFEVariant<decltype(implementedFiniteElements())>::type;

  /** \brief Default constructor
   *
   * Fills the cache with all implementations of Lagrange elements for the given
   * order and element dimension
   */
  LagrangeLocalFiniteElementCache()
  {
    cache_.resize(LocalGeometryTypeIndex::size(dim));
    Dune::Hybrid::forEach(implementedFiniteElements(), [&](auto feImpl) {
      cache_[LocalGeometryTypeIndex::index(feImpl.first)].emplace(feImpl.second());
    });
  }

  /** \brief Copy constructor */
  LagrangeLocalFiniteElementCache(const LagrangeLocalFiniteElementCache& other) = default;

  /** \brief Get the Lagrange LocalFiniteElement for the given GeometryType object
   *
   * \throws Dune::RangeError if gt has incorrect dimension
   * \throws Dune::NotImplemented if the dimension is correct, but still the cache doesn't hold
   *   a Lagrange LocalFiniteElement implementation for the requested GeometryType.
   */
  const auto& get(const GeometryType& gt) const
  {
    if (gt.dim() != dim)
      DUNE_THROW(Dune::RangeError, "You cannot get a " << gt.dim() << "-dimensional local finite element from a cache for " << dim << "-dimensional elements!");

    if (not(cache_[LocalGeometryTypeIndex::index(gt)]))
      DUNE_THROW(Dune::NotImplemented,"There is no LocalFiniteElement for a " << gt << " in the cache.");
    return *cache_[LocalGeometryTypeIndex::index(gt)];
  }

private:
  std::vector<Std::optional<FiniteElementType>> cache_;
};



} // namespace Dune




#endif // DUNE_LOCALFUNCTIONS_LAGRANGE_LAGRANGELFECACHE_HH
