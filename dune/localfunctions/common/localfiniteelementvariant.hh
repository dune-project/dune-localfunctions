// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_COMMON_LOCALFINITEELEMENTVARIANT_HH
#define DUNE_LOCALFUNCTIONS_COMMON_LOCALFINITEELEMENTVARIANT_HH

#include <cstddef>
#include <type_traits>
#include <variant>

#include <dune/common/typeutilities.hh>
#include <dune/common/std/type_traits.hh>
#include <dune/common/overloadset.hh>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>
#include <dune/localfunctions/common/localbasis.hh>
#include <dune/localfunctions/common/localkey.hh>


namespace Dune {

namespace Impl {

  // Helper for visiting a variant containing monostate.
  // Since a generic lambda will in most cases not compile
  // for monostate, we add special empty overloads for monostate.
  // Hence visitIf will simply do nothing in the case of a
  // monostate value.
  template<class Visitor, class Variant>
  void visitIf(Visitor&& visitor, Variant&& variant)
  {
    auto visitorWithFallback = overload([&](std::monostate&) {},  [&](const std::monostate&) {}, visitor);
    std::visit(visitorWithFallback, variant);
  }

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

    template<class Implementation>
    LocalBasisVariant(const Implementation& impl) :
      impl_(&impl),
      size_(impl.size()),
      order_(impl.order())
    {}

    LocalBasisVariant() = default;
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
      Impl::visitIf([&](const auto* impl) { impl->evaluateFunction(x, out); }, impl_);
    }

    /**
     * \brief Evaluate Jacobian of all shape functions
     */
    inline void evaluateJacobian(
        const typename Traits::DomainType& x,
        std::vector<typename Traits::JacobianType>& out) const
    {
      Impl::visitIf([&](const auto* impl) { impl->evaluateJacobian(x, out); }, impl_);
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
      Impl::visitIf([&](const auto* impl) { impl->partial(order, x, out); }, impl_);
    }

  private:
    std::variant<std::monostate, const Implementations*...> impl_;
    std::size_t size_;
    std::size_t order_;
  };


  template<class... Implementations>
  class LocalCoefficientsVariant
  {
  public:

    template<class Implementation>
    LocalCoefficientsVariant(const Implementation& impl) :
      impl_(&impl),
      size_(impl.size())
    {}

    LocalCoefficientsVariant() = default;
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

    const Dune::LocalKey& localKey (std::size_t i) const
    {
      // We can't use visitIf since we have to return something
      // even for a monostate value. Since the return type is
      // an l-value reference, we use a default constructed
      // dummy LocalKey value.
      static const Dune::LocalKey dummyLocalKey;
      return std::visit(overload(
          [&](const std::monostate&) -> decltype(auto) { return (dummyLocalKey);},
          [&](const auto* impl) -> decltype(auto) { return impl->localKey(i); }), impl_);
    }

  private:
    std::variant<std::monostate, const Implementations*...> impl_;
    std::size_t size_;
  };


  template<class... Implementations>
  class LocalInterpolationVariant
  {
  public:

    template<class Implementation>
    LocalInterpolationVariant(const Implementation& impl) :
      impl_(&impl)
    {}

    LocalInterpolationVariant() = default;
    LocalInterpolationVariant(const LocalInterpolationVariant& other) = default;
    LocalInterpolationVariant(LocalInterpolationVariant&& other) = default;
    LocalInterpolationVariant& operator=(const LocalInterpolationVariant& other) = default;
    LocalInterpolationVariant& operator=(LocalInterpolationVariant&& other) = default;

    template<typename F, typename C>
    void interpolate (const F& ff, std::vector<C>& out) const
    {
      Impl::visitIf([&](const auto* impl) { impl->interpolate(ff, out); }, impl_);
    }

  private:
    std::variant<std::monostate, const Implementations*...> impl_;
  };

} // namespace Impl


  /**
   * \brief Type erasure class for wrapping LocalFiniteElement classes
   *
   * This is a type erasure wrapper class for types implementing the
   * LocalFiniteElement interface. The types of the LocalFiniteElement
   * implementations that this class can hold have to be provided as
   * template parameter.
   *
   * The implementation is based on std::variant.
   * Notice that this prepends std::monostate to the Implementations
   * list for the internally stored std::variant such that
   * LocalFiniteElementVariant can be empty and is default-constructible.
   * As a consequence providing std::monostate manually to
   * LocalFiniteElementVariant is neither necessary nor allowed.
   * Access to the stored implementation is internally implemented
   * using std::visit(). To avoid multiple trivial std::visit()
   * calls, the results of size(), order(), and type() are cached
   * on creation and assignment.
   *
   * In empty state accessing any method beyond operator bool(), variant(),
   * or assignment leads to undefined behavior.
   *
   * The LocalBasisTraits are extracted from the implementation
   * provided as first template parameter. The other implementations are
   * required to be compatible with this one.
   *
   * \tparam Implementations List of supported LocalFiniteElement implementations
   */
  template<class... Implementations>
  class LocalFiniteElementVariant
  {

    // In each LocalFooVariant we store a std::variant<std::monostate, const FooImpl*...>, i.e. a std::variant
    // with the pointer to the Foo implementation unless LocalFiniteElementVariant stores a monostate. In this
    // case each LocalFooVariant also stores a monostate (and not a monostate*).
    using LocalBasis = Impl::LocalBasisVariant<typename Implementations::Traits::LocalBasisType...>;
    using LocalCoefficients = Impl::LocalCoefficientsVariant<typename Implementations::Traits::LocalCoefficientsType...>;
    using LocalInterpolation = Impl::LocalInterpolationVariant<typename Implementations::Traits::LocalInterpolationType...>;

    // Update members after changing impl_
    void updateMembers()
    {
      std::visit(overload(
          [&](std::monostate&) {
            localBasis_ = LocalBasis();
            localCoefficients_ = LocalCoefficients();
            localInterpolation_ = LocalInterpolation();
            size_ = 0;
            geometryType_ = GeometryType{};
          }, [&](auto&& impl) {
            localBasis_ = LocalBasis(impl.localBasis());
            localCoefficients_ = LocalCoefficients(impl.localCoefficients());
            localInterpolation_ = LocalInterpolation(impl.localInterpolation());
            size_ = impl.size();
            geometryType_ = impl.type();
          }), impl_);
    }

  public:

    /**
     * \brief Export LocalFiniteElementTraits
     */
    using Traits = typename Dune::LocalFiniteElementTraits<LocalBasis, LocalCoefficients, LocalInterpolation>;

    /**
     * \brief Construct empty LocalFiniteElementVariant
     */
    LocalFiniteElementVariant() = default;

    /**
     * \brief Construct empty LocalFiniteElementVariant
     */
    LocalFiniteElementVariant(const std::monostate& monostate)
    {}

    /**
     * \brief Construct LocalFiniteElementVariant
     *
     * The created LocalFiniteElementVariant will store a
     * copy of the provided implementation.
     */
    template<class Implementation,
      std::enable_if_t<std::disjunction<std::is_same<std::decay_t<Implementation>, Implementations>...>::value, int> = 0>
    LocalFiniteElementVariant(Implementation&& impl) :
      impl_(std::forward<Implementation>(impl))
    {
      updateMembers();
    }

    /**
     * \brief Copy constructor
     */
    LocalFiniteElementVariant(const LocalFiniteElementVariant& other) :
      impl_(other.impl_)
    {
      updateMembers();
    }

    /**
     * \brief Move constructor
     */
    LocalFiniteElementVariant(LocalFiniteElementVariant&& other) :
      impl_(std::move(other.impl_))
    {
      updateMembers();
    }

    /**
     * \brief Copy assignment
     */
    LocalFiniteElementVariant& operator=(const LocalFiniteElementVariant& other)
    {
      impl_ = other.impl_;
      updateMembers();
      return *this;
    }

    /**
     * \brief Move assignment
     */
    LocalFiniteElementVariant& operator=(LocalFiniteElementVariant&& other)
    {
      impl_ = std::move(other.impl_);
      updateMembers();
      return *this;
    }

    /**
     * \brief Assignment from implementation
     */
    template<class Implementation,
      std::enable_if_t<std::disjunction<std::is_same<std::decay_t<Implementation>, Implementations>...>::value, int> = 0>
    LocalFiniteElementVariant& operator=(Implementation&& impl)
    {
      impl_ = std::forward<Implementation>(impl);
      updateMembers();
      return *this;
    }


    /**
     * \brief Provide access to LocalBasis implementation of this LocalFiniteElement
     */
    const typename Traits::LocalBasisType& localBasis() const
    {
      return localBasis_;
    }

    /**
     * \brief Provide access to LocalCoefficients implementation of this LocalFiniteElement
     */
    const typename Traits::LocalCoefficientsType& localCoefficients() const
    {
      return localCoefficients_;
    }

    /**
     * \brief Provide access to LocalInterpolation implementation of this LocalFiniteElement
     */
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

    /**
     * \brief Number of shape functions
     */
    constexpr GeometryType type() const
    {
      return geometryType_;
    }

    /**
     * \brief Provide access to underlying std::variant
     *
     * This allows to use std::visit on a higher level
     * which allows to avoid the indirection of the
     * std::variant - polymorphism inside the visitor code.
     * Notice that the provided std::variant contains
     * std::monostate in its type list. Hence any
     * visitor used to access the variant has to be
     * std::monostate-aware.
     */
    const auto& variant() const
    {
      return impl_;
    }

    /**
     * \brief Check if LocalFiniteElementVariant stores an implementation
     *
     * This returns true iff variant() does not store a monostate.
     */
    operator bool () const
    {
      return not(std::holds_alternative<std::monostate>(variant()));
    }

  private:
    std::variant<std::monostate, Implementations...> impl_;
    std::size_t size_;
    GeometryType geometryType_;
    LocalBasis localBasis_;
    LocalCoefficients localCoefficients_;
    LocalInterpolation localInterpolation_;
  };

} // end namespace Dune

#endif // DUNE_LOCALFUNCTIONS_COMMON_LOCALFINITEELEMENTVARIANT_HH
