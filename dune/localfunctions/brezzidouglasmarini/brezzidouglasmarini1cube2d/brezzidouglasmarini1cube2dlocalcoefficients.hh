// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI1_CUBE2D_LOCALCOEFFICIENTS_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI1_CUBE2D_LOCALCOEFFICIENTS_HH

#include <cstddef>
#include <vector>

#include "../../common/localkey.hh"

namespace Dune
{

  /**
   * \ingroup LocalLayoutImplementation
   * \brief Layout map for Brezzi-Douglas-Marini-1 elements on quadrilaterals
   *
   * \nosubgrouping
   * \implements Dune::LocalCoefficientsVirtualImp
   */
  class BDM1Cube2DLocalCoefficients
  {

  public:
    //! \brief Standard constructor
    BDM1Cube2DLocalCoefficients () : li(8)
    {
      for (std::size_t i = 0; i < 4; ++i)
      {
        li[2*i] = LocalKey(i,1,0);
        li[2*i + 1] = LocalKey(i,1,1);
      }
    }

    //! \brief number of coefficients
    std::size_t size () const
    {
      return 8;
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
#endif // DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI1_CUBE2D_LOCALCOEFFICIENTS_HH
