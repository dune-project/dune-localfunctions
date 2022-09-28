// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS1_CUBE3D_LOCALCOEFFICIENTS_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS1_CUBE3D_LOCALCOEFFICIENTS_HH

#include <cstddef>
#include <vector>

#include "../../common/localkey.hh"

namespace Dune
{

  /**
   * \ingroup LocalLayoutImplementation
   * \brief Layout map for Raviart-Thomas-1 elements on quadrilaterals.
   *
   * \nosubgrouping
   * \implements Dune::LocalCoefficientsVirtualImp
   */
  class RT1Cube3DLocalCoefficients
  {

  public:
    //! \brief Standard constructor
    RT1Cube3DLocalCoefficients () : li(36)
    {
      for (std::size_t i = 0; i < 6; i++)
      {
        li[i] = LocalKey(i,1,0);
        li[i + 6]  = LocalKey(i,1,1);
        li[i + 12] = LocalKey(i,1,2);
        li[i + 18] = LocalKey(i,1,3);
      }

      for (std::size_t i = 0; i < 12; i++)
      {
        li[i + 24] = LocalKey(0,0,i);
      }
    }

    //! \brief number of coefficients
    std::size_t size () const
    {
      return 36;
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
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS1_CUBE3D_LOCALCOEFFICIENTS_HH
