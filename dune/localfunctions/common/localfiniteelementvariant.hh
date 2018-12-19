// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_FUNCTIONS_COMMON_LOCALFINITEELEMENTVARIANT_HH
#define DUNE_FUNCTIONS_COMMON_LOCALFINITEELEMENTVARIANT_HH

#include <cstddef>
#include <type_traits>
#include <variant>

#include <dune/common/typeutilities.hh>
#include <dune/common/std/type_traits.hh>
//#include <dune/common/std/variant.hh>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localkey.hh>


namespace Dune {

  template<class... Implementations>
  class LocalBasisVariant
  {

    template<class I0, class... II>
    struct FirstType
    { using type = I0; };

    using FirstImpTraits = typename FirstType<Implementations...>::type::Traits;

  public:

    // We do not simply copy Implementation::LocalBasisTraits because this
    // may be implementation specific. To stay clean, we simply put all its
    // data into the default LocalBasisTraits.
    using Traits = typename Dune::LocalBasisTraits<
      typename FirstImpTraits::DomainFieldType,
      FirstImpTraits::dimDomain,
      typename FirstImpTraits::DomainType,
      typename FirstImpTraits::RangeFieldType,
      FirstImpTraits::dimRange,
      typename FirstImpTraits::RangeType,
      typename FirstImpTraits::JacobianType>;

    template<class FEVariant>
    LocalBasisVariant(const FEVariant& feVariant) :
      impl_(std::visit([&](const auto& fe) { return std::variant<const Implementations*...>(&fe.localBasis()); }, feVariant)),
      size_(std::visit([&](const auto* impl) { return impl->size(); }, impl_)),
      order_(std::visit([&](const auto* impl) { return impl->order(); }, impl_))
    {}

    LocalBasisVariant() = delete;
    LocalBasisVariant(const LocalBasisVariant& other) = default;
    LocalBasisVariant(LocalBasisVariant&& other) = default;
    LocalBasisVariant& operator=(const LocalBasisVariant& other) = default;
    LocalBasisVariant& operator=(LocalBasisVariant&& other) = default;

    /**
     * \brief Number of shape functions
     */
    unsigned int size() const
    {
      return size_;
    }

    /**
     * \brief Polynomial order of the shape functions
     */
    unsigned int order() const
    {
      return order_;
    }

    /**
     * \brief Evaluate all shape functions
     */
    inline void evaluateFunction(
        const typename Traits::DomainType& x,
        std::vector<typename Traits::RangeType>& out) const
    {
      std::visit([&](const auto* impl) { impl->evaluateFunction(x, out); }, impl_);
    }

    /**
     * \brief Evaluate Jacobian of all shape functions
     */
    inline void evaluateJacobian(
        const typename Traits::DomainType& x,
        std::vector<typename Traits::JacobianType>& out) const
    {
      std::visit([&](const auto* impl) { impl->evaluateJacobian(x, out); }, impl_);
    }

    /**
     * \brief Evaluate partial derivatives of any order of all shape functions
     *
     * \param order Order of the partial derivatives, in the classic multi-index notation
     * \param x Position where to evaluate the derivatives
     * \param[out] out Return value: the desired partial derivatives
     */
    void partial(
        const std::array<unsigned int,Traits::dimDomain>& order,
        const typename Traits::DomainType& x,
        std::vector<typename Traits::RangeType>& out) const
    {
      std::visit([&](const auto* impl) { impl->partial(order, x, out); }, impl_);
    }

  private:
    std::variant<const Implementations*...> impl_;
    std::size_t size_;
    std::size_t order_;
  };


  template<class... Implementations>
  class LocalCoefficientsVariant
  {
  public:

    template<class FEVariant>
    LocalCoefficientsVariant(const FEVariant& feVariant) :
      impl_(std::visit([&](const auto& fe) { return std::variant<const Implementations*...>(&fe.localCoefficients()); }, feVariant)),
      size_(std::visit([&](const auto* impl) { return impl->size(); }, impl_))
    {}

    LocalCoefficientsVariant() = delete;
    LocalCoefficientsVariant(const LocalCoefficientsVariant& other) = default;
    LocalCoefficientsVariant(LocalCoefficientsVariant&& other) = default;
    LocalCoefficientsVariant& operator=(const LocalCoefficientsVariant& other) = default;
    LocalCoefficientsVariant& operator=(LocalCoefficientsVariant&& other) = default;

    /**
     * \brief Number of shape functions
     */
    unsigned int size() const
    {
      return size_;
    }

    Dune::LocalKey localKey (std::size_t i) const
    {
      return std::visit([&](const auto* impl) { return impl->localKey(i); }, impl_);
    }

  private:
    std::variant<const Implementations*...> impl_;
    std::size_t size_;
  };


  template<class... Implementations>
  class LocalInterpolationVariant
  {
  public:

    template<class FEVariant>
    LocalInterpolationVariant(const FEVariant& feVariant) :
      impl_(std::visit([&](const auto& fe) { return std::variant<const Implementations*...>(&fe.localInterpolation()); }, feVariant))
    {}

    LocalInterpolationVariant() = delete;
    LocalInterpolationVariant(const LocalInterpolationVariant& other) = default;
    LocalInterpolationVariant(LocalInterpolationVariant&& other) = default;
    LocalInterpolationVariant& operator=(const LocalInterpolationVariant& other) = default;
    LocalInterpolationVariant& operator=(LocalInterpolationVariant&& other) = default;

    template<typename F, typename C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      std::visit([&](const auto* impl) { impl->interpolate(ff, out); }, impl_);
    }

  private:
    std::variant<const Implementations*...> impl_;
  };



  template<class... Implementations>
  class LocalFiniteElementVariant
  {
    // In each LocalFooVariant we store a std::variant<const FooImpl*...>, i.e. a std::variant
    // with the pointer to the Foo implementation
    using LocalBasis = LocalBasisVariant<typename Implementations::Traits::LocalBasisType...>;
    using LocalCoefficients = LocalCoefficientsVariant<typename Implementations::Traits::LocalCoefficientsType...>;
    using LocalInterpolation = LocalInterpolationVariant<typename Implementations::Traits::LocalInterpolationType...>;

  public:

    using Traits = typename Dune::LocalFiniteElementTraits<LocalBasis, LocalCoefficients, LocalInterpolation>;

    LocalFiniteElementVariant() :
      impl_(),
      size_(std::visit([&](const auto& fe) { return fe.size(); }, impl_)),
      localBasis_(impl_),
      localCoefficients_(impl_),
      localInterpolation_(impl_)
    {}

    template<class Implementation,
      std::enable_if_t<Std::disjunction<std::is_same<std::decay_t<Implementation>, Implementations>...>::value, int> = 0>
    LocalFiniteElementVariant(Implementation&& impl) :
      impl_(std::forward<Implementation>(impl)),
      size_(std::visit([&](const auto& fe) { return fe.size(); }, impl_)),
      localBasis_(impl_),
      localCoefficients_(impl_),
      localInterpolation_(impl_)
    {}

    LocalFiniteElementVariant(const LocalFiniteElementVariant& other) :
      impl_(other.impl_),
      size_(other.size_),
      localBasis_(impl_),
      localCoefficients_(impl_),
      localInterpolation_(impl_)
    {}

    LocalFiniteElementVariant(LocalFiniteElementVariant&& other) :
      impl_(std::move(other.impl_)),
      size_(other.size_),
      localBasis_(impl_),
      localCoefficients_(impl_),
      localInterpolation_(impl_)
    {}

    LocalFiniteElementVariant& operator=(const LocalFiniteElementVariant& other)
    {
      impl_ = other.impl_;
      size_ = other.size_;
      localBasis_ = LocalBasis(impl_);
      localCoefficients_ = LocalCoefficients(impl_);
      localInterpolation_ = LocalInterpolation(impl_);
      return *this;
    }

    const typename Traits::LocalBasisType& localBasis() const
    {
      return localBasis_;
    }

    const typename Traits::LocalCoefficientsType& localCoefficients() const
    {
      return localCoefficients_;
    }

    const typename Traits::LocalInterpolationType& localInterpolation() const
    {
      return localInterpolation_;
    }

    /**
     * \brief Number of shape functions
     */
    unsigned int size() const
    {
      return size_;
    }

    constexpr GeometryType type() const
    {
      return std::visit([&](const auto& fe) { return fe.type(); }, impl_);
    }

    /**
     * \brief Provide access to undelying std::variant
     *
     * This allows to use std::visit on a higher level
     * which allows to avoid the indirection of the
     * std::variant - polymorphism inside the visitor code.
     */
    const auto& variant() const
    {
      return impl_;
    }

  private:
    std::variant<Implementations...> impl_;
    std::size_t size_;
    LocalBasis localBasis_;
    LocalCoefficients localCoefficients_;
    LocalInterpolation localInterpolation_;
  };

} // end namespace Dune

#endif // DUNE_FUNCTIONS_COMMON_LOCALFINITEELEMENTVARIANT_HH
