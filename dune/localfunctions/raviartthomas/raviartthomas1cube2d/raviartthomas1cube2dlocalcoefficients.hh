// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS1_CUBE2D_LOCALCOEFFICIENTS_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS1_CUBE2D_LOCALCOEFFICIENTS_HH

#include <cstddef>
#include <vector>

#include "../../common/localkey.hh"

namespace Dune
{

  /**
   * \ingroup LocalLayoutImplementation
   * \brief Layout map for Raviart-Thomas-1 elements on quadrilaterals
   *
   * \nosubgrouping
   * \implements Dune::LocalCoefficientsVirtualImp
   */
  class RT1Cube2DLocalCoefficients
  {

  public:
    //! \brief Standard constructor
    RT1Cube2DLocalCoefficients () : li(12)
    {
      for (std::size_t i=0; i < 4; i++)
      {
        li[2*i] = LocalKey(i,1,0);
        li[2*i + 1] = LocalKey(i,1,1);
      }

      li[8]  = LocalKey(0,0,0);
      li[9]  = LocalKey(0,0,1);
      li[10] = LocalKey(0,0,2);
      li[11] = LocalKey(0,0,3);
    }

    //! number of coefficients
    std::size_t size () const
    {
      return 12;
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
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS1_CUBE2D_LOCALCOEFFICIENTS_HH
