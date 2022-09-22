// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_COMMON_LOCALBASIS_HH
#define DUNE_LOCALFUNCTIONS_COMMON_LOCALBASIS_HH

namespace Dune
{

  /**@ingroup LocalBasisInterface
         \brief Type traits for LocalBasisVirtualInterface

         A shape function is a function
         \f[ \hat\phi : \mbox{IR}^n \to \mbox{IR}^m. \f]
         This traits class holds information how the signature of this
         function is represented in C++ types.

         This is just a convenience class for supplying traits to the
         LocalBasisVirtualInterface and its implementations.

         \tparam DF Type to represent the field in the domain.
         \tparam n  Dimension of the domain.
         \tparam D  Type to represent the domain, allows random access.
         \tparam RF Type to represent the field in the range.
         \tparam m  Dimension of the range.
         \tparam R  Type to represent the range, allows random access.
         \tparam J  Type to represent the Jacobian, allows random access.

         \nosubgrouping
   */
  template<class DF, int n, class D, class RF, int m, class R, class J>
  struct LocalBasisTraits
  {
    //! \brief Export type for domain field
    typedef DF DomainFieldType;

    //! \brief dimension of the domain
    constexpr static int dimDomain = n;

    //! \brief domain type
    typedef D DomainType;

    //! \brief Export type for range field
    typedef RF RangeFieldType;

    //! \brief dimension of the range
    constexpr static int dimRange = m;

    //! \brief range type
    typedef R RangeType;

    /** \brief Type to represent derivative

            When \f$ \hat\phi : \mbox{IR}^n \to \mbox{IR}^m \f$ then JacobianType
            is an 2D-array of m x n components where entry J[i][j] contains
            the derivative  \f$\partial_j \hat\phi_i \f$.
     */
    typedef J JacobianType;
  };

}
#endif
