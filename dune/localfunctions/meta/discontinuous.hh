// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#ifndef DUNE_LOCALFUNCTIONS_META_DISCONTINUOUS_HH
#define DUNE_LOCALFUNCTIONS_META_DISCONTINUOUS_HH

#include <utility>

#include <dune/common/referencehelper.hh>
#include <dune/geometry/type.hh>
#include <dune/localfunctions/utility/dglocalcoefficients.hh>

namespace Dune {

  //! \brief Meta-finite element turning a finite-element into "discontinuous" finite-element
  //! by associating all basis functions to the element interior.
  /**
   * \ingroup LocalFunctions
   *
   * \tparam LFE  Type of the local finite-element, can be a raw type or a reference_wrapper.
   */
  template<class LFE>
  class DiscontinuousLocalFiniteElement
  {
    using LFERaw = ResolveRef_t<LFE>;

    using LB = typename LFERaw::Traits::LocalBasisType;
    using LC = typename LFERaw::Traits::LocalCoefficientsType;
    using LI = typename LFERaw::Traits::LocalInterpolationType;

  public:
    //! types of component objects
    struct Traits {
      //! type of the Basis
      using LocalBasisType = LB;
      //! type of the Coefficients
      using LocalCoefficientsType = DGLocalCoefficients;
      //! type of the Interpolation
      using LocalInterpolationType = LI;
    };

  private:
    LFE lfe_;
    typename Traits::LocalCoefficientsType lc_;

  public:
    //! Construct a finite element
    template <class LFE_>
    explicit DiscontinuousLocalFiniteElement (LFE_&& lfe)
      : lfe_(std::forward<LFE_>(lfe))
      , lc_(resolveRef(lfe_).localCoefficients().size())
    {}

    //! Extract basis of this finite element
    /**
     * The returned lvalue must have a lifetime at least as long as the finite
     * element object it was acquired from.
     */
    const typename Traits::LocalBasisType& localBasis () const
    {
      return resolveRef(lfe_).localBasis();
    }

    //! Extract coefficients of this finite element
    /**
     * The returned lvalue must have a lifetime at least as long as the finite
     * element object it was acquired from.
     */
    const typename Traits::LocalCoefficientsType& localCoefficients() const
    {
      return lc_;
    }

    //! Extract interpolation of this finite element
    /**
     * The returned lvalue must have a lifetime at least as long as the finite
     * element object it was acquired from.
     */
    const typename Traits::LocalInterpolationType& localInterpolation() const
    {
      return resolveRef(lfe_).localInterpolation();
    }

    //! Return the number of basis functions
    unsigned int size () const
    {
      return resolveRef(lfe_).size();
    }

    //! Return the geometry type the finite element can be bound to
    const GeometryType type () const
    {
      return resolveRef(lfe_).type();
    }
  };

  // deduction guide
  template <class LFE>
  DiscontinuousLocalFiniteElement (LFE lfe)
    -> DiscontinuousLocalFiniteElement<LFE>;

} // namespace Dune

#endif // DUNE_LOCALFUNCTIONS_META_DISCONTINUOUS_HH
