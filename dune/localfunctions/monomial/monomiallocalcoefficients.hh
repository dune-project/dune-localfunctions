// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_MONOMIAL_MONOMIALLOCALCOEFFICIENTS_HH
#define DUNE_LOCALFUNCTIONS_MONOMIAL_MONOMIALLOCALCOEFFICIENTS_HH

#include <cstddef>
#include <vector>

#include "../common/localkey.hh"

namespace Dune
{

  /**@ingroup LocalLayoutImplementation
         \brief Layout map for monomial finite elements

         \nosubgrouping
     \implements Dune::LocalCoefficientsVirtualImp
   */
  template<int static_size>
  class MonomialLocalCoefficients
  {
  public:
    //! \brief Standard constructor
    MonomialLocalCoefficients ()
      : index(static_size, LocalKey(0,0,0))
    {
      for(int i = 0; i < static_size; ++i)
        index[i].index(i);
    }

    //! number of coefficients
    std::size_t size () const
    {
      return static_size;
    }

    //! get i'th index
    const LocalKey& localKey (std::size_t i) const
    {
      return index[i];
    }

  private:
    std::vector<LocalKey> index;
  };

}
#endif //DUNE_LOCALFUNCTIONS_MONOMIAL_MONOMIALLOCALCOEFFICIENTS_HH
