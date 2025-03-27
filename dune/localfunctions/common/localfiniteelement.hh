// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_COMMON_TYPEERASEDLOCALFINITEELEMENT
#define DUNE_LOCALFUNCTIONS_COMMON_TYPEERASEDLOCALFINITEELEMENT

#include <memory>
#include <type_traits>
#include <utility>

#include <dune/common/typeutilities.hh>

#include <dune/geometry/type.hh>


#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>

namespace Dune
{

  /**
   * \brief Type erasure class storing a local finite element
   *
   * \tparam LBT A LocalBasisTraits class encoding the local basis signature.
   *
   * In the flavor of `std::function` this class implements a type
   * erased wrapper around local finite element implementations.
   * One can either construct the wrapper directly with the
   * local finite element that should be stored or construct
   * and empty wrapper and assign it later.
   * The type erasure is internally implemented using LocalFiniteElementVirtualInterface
   * and LocalFiniteElementVirtualImp. In contrast to these,
   * the LocalFiniteElement provides value-semantics and encapsulates
   * dynamic memory management.
   */
  template<class LBT>
  class LocalFiniteElement
  {
    using LocalBasisTraits = LBT;
    using VirtualLFE = typename Dune::LocalFiniteElementVirtualInterface<LBT>;

  public:

    using Traits = typename VirtualLFE::Traits;

    /**
     * \brief Default constructor
     *
     * The constructed wrapper is empty and can later
     * be filled with a local finite element using the
     * assignment operators.
     */
    LocalFiniteElement() = default;

    /**
     * \brief Construct from implementation
     *
     * \tparam LFEImpl LocalFiniteElement implementation type
     *
     * \param lfe LocalFiniteElement implementation
     *
     * Construct a LocalFiniteElement wrapper storing a copy
     * of the given local finite element object.
     */
    template<class LFEImpl,
      Dune::disableCopyMove<LocalFiniteElement, LFEImpl> = 0 >
    LocalFiniteElement(LFEImpl&& lfe)
      : lfe_(new LocalFiniteElementVirtualImp<std::decay_t<LFEImpl>>(std::forward<LFEImpl>(lfe)))
      , size_(lfe.size())
      , type_(lfe.type())
    {}

    /**
     * \brief Copy constructor
     *
     * \param other The LocalFiniteElement to copy from
     *
     * Depending of the state of the argument, the newly constructed
     * LocalFiniteElement will either be empty or contain a copy of the
     * local finite element stored in the LocalFiniteElement passed as
     * constructor argument.
     */
    LocalFiniteElement(const LocalFiniteElement& other)
      : lfe_(other.lfe_ ? other.lfe_->clone() : nullptr)
      , size_(other.size())
      , type_(other.type())
    {}

    /**
     * \brief Move constructor
     *
     * \param other The LocalFiniteElement to move from
     *
     * Depending of the state of the argument, the newly constructed
     * LocalFiniteElement will either be empty or contain the
     * local finite element that was originally stored in the
     * LocalFiniteElement passed as constructor argument.
     */
    LocalFiniteElement(LocalFiniteElement&& other) = default;

    /**
     * \brief Copy assignment
     *
     * \param rhs The LocalFiniteElement to assign from
     *
     * This will first remove the stored local finite element. Then,
     * depending of the state of the argument, `*this` will either
     * be set to empty state or store a copy of the local finite element
     * stored in the LocalFiniteElement passed as right hand side.
     */
    LocalFiniteElement& operator=(const LocalFiniteElement& rhs)
    {
      if (&rhs!=this)
      {
        if (rhs.lfe_)
          lfe_.reset(rhs.lfe_->clone());
        else
          lfe_.release();
        size_ = rhs.size();
        type_ = rhs.type();
      }
      return *this;
    }

    /**
     * \brief Move assignment
     *
     * \param rhs The LocalFiniteElement to assign from
     *
     * This will first remove the stored local finite element. Then,
     * depending of the state of the argument, `*this` will either
     * be set to empty state or take over the local finite element
     * stored in the LocalFiniteElement passed as right hand side.
     * In any case the latter will be empty afterwards.
     */
    LocalFiniteElement& operator=(LocalFiniteElement&& other) = default;

    /**
     * \brief Check if the LocalFiniteElement is not empty
     *
     * This method returns `false` if `*this` is empty and `true` otherwise
     */
    operator bool () const
    {
      return lfe_;
    }

    /**
     * \brief Access the LocalBasis of the stored local finite element
     *
     * Calling this method if the LocalFiniteElement is empty
     * will invoke undefined behavior.
     */
    const typename Traits::LocalBasisType& localBasis () const
    {
      return lfe_->localBasis();
    }

    /**
     * \brief Access the LocalCoefficients of the stored local finite element
     *
     * Calling this method if the LocalFiniteElement is empty
     * will invoke undefined behavior.
     */
    const typename Traits::LocalCoefficientsType& localCoefficients () const
    {
      return lfe_->localCoefficients();
    }

    /**
     * \brief Access the LocalInterpolation of the stored local finite element
     *
     * Calling this method if the LocalFiniteElement is empty
     * will invoke undefined behavior.
     */
    const typename Traits::LocalInterpolationType& localInterpolation () const
    {
      return lfe_->localInterpolation();
    }

    /**
     * \brief Get the number of basis functions of the stored local finite element
     *
     * The result of this method is unspecified if the LocalFiniteElement is empty.
     */
    unsigned int size () const
    {
      return size_;
    }

    /**
     * \brief Get the GeometryType of the stored local finite element
     *
     * The result of this method is unspecified if the LocalFiniteElement is empty.
     */
    const GeometryType& type () const
    {
      return type_;
    }

  private:
    std::unique_ptr<VirtualLFE> lfe_;
    unsigned int size_ = 0;
    Dune::GeometryType type_;
  };

}
#endif
