// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_COMMON_LOCALFINITEELEMENTVARIANTCACHE_HH
#define DUNE_LOCALFUNCTIONS_COMMON_LOCALFINITEELEMENTVARIANTCACHE_HH

#include <vector>
#include <tuple>
#include <utility>
#include <type_traits>

#include <dune/common/std/type_traits.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/typelist.hh>
#include <dune/common/hybridutilities.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

#include <dune/localfunctions/common/localfiniteelementvariant.hh>


namespace Dune {

namespace Impl {

  // This class provides the index method of LocalGeometryTypeIndex
  // but throws a Dune::RangeError if the dimension does not match.
  // This can be helpful to catch errors in a LocalFiniteElementVariantCache
  // instance based on dimension specific GeometryType indices.
  template<std::size_t dim>
  struct FixedDimLocalGeometryTypeIndex {
    inline static std::size_t index(const GeometryType &gt)
    {
      if (gt.dim() != dim)
        DUNE_THROW(Dune::RangeError, "Asking for dim=" << dim << " specific index of GeometryType with dimension " << gt.dim());
      return LocalGeometryTypeIndex::index(gt);
    }
  };

} // end namespace Impl

/** \brief A cache storing a compile time selection of local finite element implementations
 *
 * This class stores a compile time selection of LocalFiniteElement implementations.
 * The exported FiniteElementType is a LocalFiniteElementVariant<Implementations...>.
 *
 * The base class is required to implement two methods: getImplementations() and index().
 * The index(key...) method computes an index that uniquely identifies a LocalFiniteElement
 * implementation for the provided key data.
 * The getImplementations() method returns a std::tuple with one entry for each
 * LocalFiniteElement implementation. Each entry is an std::pair consisting of the
 * index of the implementation and a callable object that creates a corresponding
 * LocalFiniteElement implementation.
 *
 * The constructor forwards all arguments to the base class.
 *
 * For performence reasons all LocalFiniteElement implementations are created
 * during construction and stored in an std::vector according to the associated
 * indices. Hence densly packed indices should be prefered to avoid wasting space
 * for emplty LocalFiniteElementVariant's string no implementation.
 *
 * \tparam Base Type of the base class providing getImplementations() and index().
 */
template<class Base>
class LocalFiniteElementVariantCache : Base
{

  template<class LFEImplTuple>
  struct GenerateLFEVariant;

  template<class Index, class... LFEImpl>
  struct GenerateLFEVariant<std::tuple<std::pair<Index, LFEImpl>...>>
  {
    using type = UniqueTypes_t<LocalFiniteElementVariant, decltype(std::declval<LFEImpl>()())...>;
  };

  using Base::getImplementations;
  using Base::index;
  using Implementations = decltype(std::declval<Base>().getImplementations());

public:

  /**
   * \brief Type of exported LocalFiniteElement's
   *
   * This is a LocalFiniteElementVariant<Implementation...>
   * with Implementations... being the list of all provided
   * implementations.
   */
  using FiniteElementType = typename GenerateLFEVariant<Implementations>::type;

  /** \brief Default constructor
   *
   * Prefills the cache with all implementations.
   */
  template<class... Args>
  LocalFiniteElementVariantCache(Args&&... args) :
    Base(std::forward<Args>(args)...)
  {
    Dune::Hybrid::forEach(getImplementations(), [&,this](auto feImpl) {
      auto implIndex = feImpl.first;
      if (cache_.size() < implIndex+1)
        cache_.resize(implIndex+1);
      cache_[implIndex] = feImpl.second();
    });
  }

  /** \brief Copy constructor */
  LocalFiniteElementVariantCache(const LocalFiniteElementVariantCache& other) = default;

  /** \brief Move constructor */
  LocalFiniteElementVariantCache(LocalFiniteElementVariantCache&& other) = default;

  /** \brief Get the LocalFiniteElement for the given key data
   *
   * \throws Dune::RangeError If the cache doesn't hold a value matching the requested type.
   */
  template<class... Key>
  const auto& get(const Key&... key) const
  {
    auto implIndex = index(key...);
    if (implIndex >= cache_.size())
      DUNE_THROW(Dune::RangeError,"There is no LocalFiniteElement of the requested type.");
    if (not(cache_[implIndex]))
      DUNE_THROW(Dune::RangeError,"There is no LocalFiniteElement of the requested type.");
    return cache_[implIndex];
  }

private:
  std::vector<FiniteElementType> cache_;
};



} // namespace Dune




#endif // DUNE_LOCALFUNCTIONS_COMMON_LOCALFINITEELEMENTVARIANT_HH
