// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_RANNACHER_TUREK_LOCALFINITEELEMENT_HH
#define DUNE_RANNACHER_TUREK_LOCALFINITEELEMENT_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localfiniteelementtraits.hh>

#include "rannachertureklocalbasis.hh"
#include "rannachertureklocalcoefficients.hh"
#include "rannachertureklocalinterpolation.hh"

namespace Dune
{

  /**
   * \brief Rannacher-Turek shape functions.
   *
   * \ingroup RannacherTurek
   *
   * \tparam D type to represent the field in the domain.
   * \tparam R type to represent the field in the range.
   * \tparam d domain dimension
   */
  template< class D, class R, unsigned int d >
  struct RannacherTurekLocalFiniteElement
  {
    //! \brief export traits class
    typedef LocalFiniteElementTraits< RannacherTurekLocalBasis< D, R, d >,
        RannacherTurekLocalCoefficients< d >,
        RannacherTurekLocalInterpolation< D, R, d >
        > Traits;

    //! \brief return local basis
    const typename Traits::LocalBasisType &localBasis () const
    {
      return localBasis_;
    }

    //! \brief return local coefficients
    const typename Traits::LocalCoefficientsType &localCoefficients () const
    {
      return localCoefficients_;
    }

    //! \brief return local interpolation
    const typename Traits::LocalInterpolationType &localInterpolation () const
    {
      return localInterpolation_;
    }

    /** \brief Number of shape functions in this finite element */
    unsigned int size () const
    {
      return localBasis_.size();
    }

    //! \brief return geometry type
    GeometryType type () const
    {
      return GeometryTypes::cube(d);
    }

  private:
    typename Traits::LocalBasisType localBasis_;
    typename Traits::LocalCoefficientsType localCoefficients_;
    typename Traits::LocalInterpolationType localInterpolation_;
  };

} // namespace Dune

#endif // #ifndef DUNE_RANNACHER_TUREK_LOCALFINITEELEMENT_HH
