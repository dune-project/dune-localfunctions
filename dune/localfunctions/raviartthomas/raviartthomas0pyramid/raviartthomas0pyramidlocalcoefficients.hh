// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_PYRAMID_LOCALCOEFFICIENTS_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_PYRAMID_LOCALCOEFFICIENTS_HH

#include <cstddef>
#include <vector>

#include "../../common/localkey.hh"

namespace Dune
{

  /**
   * \brief Layout map for Raviart-Thomas-1 elements on pyramids.
   *
   * \nosubgrouping
   * \ingroup RaviartThomasImpl
   * \implements Dune::LocalCoefficientsVirtualImp
   */
  class RT0PyramidLocalCoefficients
  {

  public:
    //! \brief Standard constructor
    RT0PyramidLocalCoefficients () : li(size())
    {
      for(std::size_t i=0; i< size(); i++)
        li[i] = LocalKey(i,1,0);
    }

    //! \brief number of coefficients
    std::size_t size () const
    {
      return 5;
    }

    //! \brief get i'th index
    const LocalKey& localKey (std::size_t i) const
    {
      return li[i];
    }

  private:
    std::vector<LocalKey> li;
  };
}
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS0_PYRAMID_LOCALCOEFFICIENTS_HH
