// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_RANNACHER_TUREK_LOCALCOEFFICIENTS_HH
#define DUNE_RANNACHER_TUREK_LOCALCOEFFICIENTS_HH

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>

#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{
  /**@ingroup LocalLayoutImplementation
     \brief layout for Rannacher-Turek elements

     \tparam d Domain dimension

     \nosubgrouping
   */
  template< unsigned int d >
  struct RannacherTurekLocalCoefficients
  {
    RannacherTurekLocalCoefficients ()
    {
      for( std::size_t i = 0; i < 2*d; ++i )
        localKeys_[ i ] = LocalKey( i, 1, 0 );
    }

    RannacherTurekLocalCoefficients ( const RannacherTurekLocalCoefficients &other )
    {
      (*this) = other;
    }

    RannacherTurekLocalCoefficients &operator= ( const RannacherTurekLocalCoefficients &other )
    {
      std::copy( other.localKeys_.begin(), other.localKeys_.end(), localKeys_.begin() );
      return *this;
    }

    //! number of coefficients
    std::size_t size () const
    {
      return 2*d;
    }

    //! map index i to local key
    const LocalKey &localKey ( std::size_t i ) const
    {
      assert( i < 2*d );
      return localKeys_[ i ];
    }

  private:
    std::array< LocalKey, 2*d > localKeys_;
  };

} // namespace Dune

#endif // #ifndef DUNE_RANNACHER_TUREK_LOCALCOEFFICIENTS_HH
