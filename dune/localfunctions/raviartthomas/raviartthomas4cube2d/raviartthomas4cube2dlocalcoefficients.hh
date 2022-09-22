// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS4_CUBE2D_LOCALCOEFFICIENTS_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS4_CUBE2D_LOCALCOEFFICIENTS_HH

#include <cstddef>
#include <vector>

#include "../../common/localkey.hh"

namespace Dune
{
  /**
   * \ingroup LocalLayoutImplementation
   * \brief Layout map for Raviart-Thomas-4 elements on quadrilaterals
   *
   * \nosubgrouping
   * \implements Dune::LocalCoefficientsVirtualImp
   */
  class RT4Cube2DLocalCoefficients
  {

  public:
    //! \brief Standard constructor
    RT4Cube2DLocalCoefficients () : li(60)
    {
      for (std::size_t i = 0; i < 4; i++)
      {
        li[5*i] = LocalKey(i,1,0);
        li[5*i + 1] = LocalKey(i,1,1);
        li[5*i + 2] = LocalKey(i,1,2);
        li[5*i + 3] = LocalKey(i,1,3);
        li[5*i + 4] = LocalKey(i,1,4);
      }

      for (std::size_t i=0; i<40; i++)
      {
        li[20 + i] = LocalKey(0,0,i);
      }
    }

    //! \brief number of coefficients
    std::size_t size () const
    {
      return 60;
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
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS4_CUBE2D_LOCALCOEFFICIENTS_HH
