// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_P0LOCALCOEFFICIENTS_HH
#define DUNE_P0LOCALCOEFFICIENTS_HH

#include <cstddef>
#include <iostream>
#include <vector>

#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{

  /**@ingroup LocalLayoutImplementation
         \brief Layout map for P0 elements

         \nosubgrouping
     \implements Dune::LocalCoefficientsVirtualImp
   */
  class P0LocalCoefficients
  {
  public:
    //! \brief Standard constructor
    P0LocalCoefficients () : index(0,0,0)
    {}

    //! number of coefficients
    std::size_t size () const
    {
      return 1;
    }

    //! get i'th index
    const LocalKey& localKey ([[maybe_unused]] std::size_t i) const
    {
      return index;
    }

  private:
    LocalKey index;
  };

}
#endif
