// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS_RAVIARTTHOMASLFECACHE_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS_RAVIARTTHOMASLFECACHE_HH

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

#include <dune/localfunctions/common/localfiniteelementvariant.hh>
#include <dune/localfunctions/common/uniquetypes.hh>
#include <dune/localfunctions/raviartthomas.hh>


namespace Dune {

namespace Impl {

  // Provide implemented Raviart-Thomas local finite elements

  template<class D, class R, std::size_t dim, std::size_t order>
  struct ImplementedRaviartThomasLocalFiniteElements
  {};

  template<class D, class R>
  struct ImplementedRaviartThomasLocalFiniteElements<D,R,2,0>
  {
    static auto get()
    {
      return std::make_tuple(
        std::make_pair(GeometryTypes::triangle, []() { return RT02DLocalFiniteElement<D,R>(); }),
        std::make_pair(GeometryTypes::quadrilateral, []() { return RT0Cube2DLocalFiniteElement<D,R>(); })
      );
    }
  };

  template<class D, class R>
  struct ImplementedRaviartThomasLocalFiniteElements<D,R,2,1>
  {
    static auto get()
    {
      return std::make_tuple(
        std::make_pair(GeometryTypes::triangle, []() { return RT12DLocalFiniteElement<D,R>(); }),
        std::make_pair(GeometryTypes::quadrilateral, []() { return RT1Cube2DLocalFiniteElement<D,R>(); })
      );
    }
  };

  template<class D, class R>
  struct ImplementedRaviartThomasLocalFiniteElements<D,R,2,2>
  {
    static auto get()
    {
      return std::make_tuple(
        std::make_pair(GeometryTypes::quadrilateral, []() { return RT2Cube2DLocalFiniteElement<D,R>(); })
      );
    }
  };

  template<class D, class R>
  struct ImplementedRaviartThomasLocalFiniteElements<D,R,3,0>
  {
    static auto get()
    {
      return std::make_tuple(
        std::make_pair(GeometryTypes::hexahedron, []() { return RT0Cube3DLocalFiniteElement<D,R>(); })
      );
    }
  };

  template<class D, class R>
  struct ImplementedRaviartThomasLocalFiniteElements<D,R,3,1>
  {
    static auto get()
    {
      return std::make_tuple(
        std::make_pair(GeometryTypes::hexahedron, []() { RT1Cube3DLocalFiniteElement<D,R>(); })
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
 */
template<class D, class R, std::size_t dim, std::size_t order>
class RaviartThomasLocalFiniteElementCache
{
  template<class LFEImplTuple>
  struct GenerateLFEVariant;

  template<class... LFEImpl>
  struct GenerateLFEVariant<std::tuple<std::pair<GeometryType, LFEImpl>...>>
  {
    using type = UniqueTypes_t<LocalFiniteElementVariant, decltype(std::declval<LFEImpl>()())...>;
  };

  static auto implementedLocalFiniteElements()
  {
    return Impl::ImplementedRaviartThomasLocalFiniteElements<D,R,dim, order>::get();
  }

public:

  using FiniteElementType = typename GenerateLFEVariant<decltype(implementedLocalFiniteElements())>::type;

  /** \brief Default constructor
   *
   * Fills the cache with all implementations of Raviart-Thomas elements for the given
   * order and element dimension
   */
  RaviartThomasLocalFiniteElementCache()
  {
    cache_.resize(LocalGeometryTypeIndex::size(dim));
    Dune::Hybrid::forEach(implementedLocalFiniteElements(), [&](auto feImpl) {
      cache_[LocalGeometryTypeIndex::index(feImpl.first)].emplace(feImpl.second());
    });
  }

  /** \brief Copy constructor */
  RaviartThomasLocalFiniteElementCache(const RaviartThomasLocalFiniteElementCache& other) = default;

  /** \brief Get the RaviartThomas LocalFiniteElement for the given GeometryType object
   *
   * \throws Dune::RangeError if gt has incorrect dimension
   * \throws Dune::NotImplemented if the dimension is correct, but still the cache doesn't hold
   *   a Raviart-Thomas LocalFiniteElement implementation for the requested GeometryType.
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

#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS_RAVIARTTHOMASLFECACHE_HH
