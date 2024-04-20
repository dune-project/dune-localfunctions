// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

#ifndef DUNE_LOCALFUNCTIONS_META_DISCONTINUOUS_HH
#define DUNE_LOCALFUNCTIONS_META_DISCONTINUOUS_HH

#include <memory>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/utility/dglocalcoefficients.hh>

namespace Dune {

  //! \brief Meta-finite element turning a finite-element into "discontinuous" finite-element
  //! by associating all basis functions to the element interior.
  /**
   * \ingroup LocalFunctions
   *
   * \tparam LFE  Type of the local finite-element
   */
  template<class LFE>
  class DiscontinuousLocalFiniteElement
  {
    using LB = typename LFE::Traits::LocalBasisType;
    using LC = typename LFE::Traits::LocalCoefficientsType;
    using LI = typename LFE::Traits::LocalInterpolationType;

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
    std::shared_ptr<const LFE> lfe_;
    typename Traits::LocalCoefficientsType lc_;

  public:
    //! Construct a finite element
    DiscontinuousLocalFiniteElement (const LFE& lfe)
      : DiscontinuousLocalFiniteElement(std::make_shared<const LFE>(lfe))
    {}

    //! Construct a finite element
    DiscontinuousLocalFiniteElement (std::shared_ptr<const LFE> lfe)
      : lfe_(std::move(lfe))
      , lc_(lfe_->localCoefficients().size())
    {}

    //! Extract basis of this finite element
    /**
     * The returned lvalue must have a lifetime at least as long as the finite
     * element object it was acquired from.
     */
    const typename Traits::LocalBasisType& localBasis () const
    {
      return lfe_->localBasis();
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
      return lfe_->localInterpolation();
    }

    //! Return the number of basis functions
    unsigned int size () const
    {
      return lfe_->size();
    }

    //! Return the prismatic-product of the lfe's geometry types
    const GeometryType type () const
    {
      return lfe_->type();
    }
  };

} // namespace Dune

#endif // DUNE_LOCALFUNCTIONS_META_DISCONTINUOUS_HH
