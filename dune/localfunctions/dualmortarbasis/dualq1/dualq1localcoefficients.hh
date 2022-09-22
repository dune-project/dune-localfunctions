// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_DUAL_Q1_LOCALCOEFFICIENTS_HH
#define DUNE_DUAL_Q1_LOCALCOEFFICIENTS_HH

#include <cstddef>
#include <iostream>
#include <vector>

#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{

  /**@ingroup LocalLayoutImplementation
         \brief Layout map for dual Q1 elements

         \nosubgrouping
     \implements Dune::LocalCoefficientsVirtualImp
   */
  template <int dim>
  class DualQ1LocalCoefficients
  {
  public:
    //! \brief Standard constructor
    DualQ1LocalCoefficients () : li(1<<dim)
    {
      for (std::size_t i=0; i<(1<<dim); i++)
        li[i] = LocalKey(i,dim,0);
    }

    //! number of coefficients
    std::size_t size () const
    {
      return 1<<dim;
    }

    //! get i'th index
    const LocalKey& localKey (std::size_t i) const
    {
      return li[i];
    }

  private:
    std::vector<LocalKey> li;
  };

}

#endif
