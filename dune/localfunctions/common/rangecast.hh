// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_COMMON_RANGECAST_HH
#define DUNE_LOCALFUNCTIONS_COMMON_RANGECAST_HH

#include <dune/common/typetraits.hh>
#include <dune/common/fvector.hh>



namespace Dune {

namespace Impl {

  /**
   * \brief Cast LocalBasis range type
   *
   * This overload is selected if the argument
   * is a number type and forwards the argument.
   */
  template<class X>
  requires (Dune::IsNumber<X>::value)
  const X& rangeCast(const X& x)
  {
    return x;
  }

  /**
   * \brief Cast LocalBasis range type
   *
   * This overload is selected for FieldVector<X,1>. This type
   * can depending on the application either be used to represent
   * a number or a single-entry vector. In any case casting to
   * the entry is fine, because we can assign a FieldVector<X,1>
   * from X.
   */
  template<class X>
  const X& rangeCast(const Dune::FieldVector<X, 1>& x)
  {
    return x[0];
  }

  /**
   * \brief Cast LocalBasis range type
   *
   *
   * This overload is selected for FieldVector<X,n> with n>1.
   * This type is treated as a vector and hence forwarded
   * directly.
   */
  template<class X, int n>
  requires(n>1)
  const X& rangeCast(const Dune::FieldVector<X, n>& x)
  {
    return x;
  }

} // end namespace Impl

} // end namespace Dune



#endif // DUNE_LOCALFUNCTIONS_COMMON_RANGECAST_HH
