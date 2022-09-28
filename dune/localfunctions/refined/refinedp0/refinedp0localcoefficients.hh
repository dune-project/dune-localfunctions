// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_REFINED_P0_LOCALCOEFFICIENTS_HH
#define DUNE_REFINED_P0_LOCALCOEFFICIENTS_HH

#include <cstddef>
#include <iostream>
#include <vector>

#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{

  /**@ingroup LocalLayoutImplementation
     \brief Layout map for RefinedP0 elements

     \nosubgrouping
     \implements Dune::LocalCoefficientsVirtualImp
   */
  template<unsigned int k>
  class RefinedP0LocalCoefficients
  {
    // 2 to the k-th power
    constexpr static int N = 1<<k;

  public:
    RefinedP0LocalCoefficients () :
      localKeys_(N)
    {
      // All functions are associated to the element
      for (int i = 0; i < N; ++i)
        localKeys_[i] = LocalKey(0,0,i);
    }

    //! number of coefficients
    std::size_t size () const
    {
      return N;
    }

    //! get i'th index
    const LocalKey& localKey (std::size_t i) const
    {
      return localKeys_[i];
    }

  private:
    std::vector<LocalKey> localKeys_;

  };

}

#endif
