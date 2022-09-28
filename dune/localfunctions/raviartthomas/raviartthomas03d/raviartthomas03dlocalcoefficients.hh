// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS_RAVIARTTHOMAS03D_RAVIARTTHOMAS03DLOCALCOEFFICIENTS_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS_RAVIARTTHOMAS03D_RAVIARTTHOMAS03DLOCALCOEFFICIENTS_HH

#include <cstddef>
#include <iostream>
#include <vector>

#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{

  /**@ingroup LocalLayoutImplementation
         \brief Layout map for RT0 elements

         \nosubgrouping
     \implements Dune::LocalCoefficientsVirtualImp
   */
  class RT03DLocalCoefficients
  {
  public:
    //! \brief Standard constructor
    RT03DLocalCoefficients () : li(4)
    {
      for (std::size_t i=0; i<4; i++)
        li[i] = LocalKey(i,1,0);
    }

    //! number of coefficients
    std::size_t size () const
    {
      return 4;
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
