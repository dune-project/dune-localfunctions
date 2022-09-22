// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
#ifndef DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI2_CUBE2D_LOCALCOEFFICIENTS_HH
#define DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI2_CUBE2D_LOCALCOEFFICIENTS_HH

#include <cstddef>
#include <vector>

#include "../../common/localkey.hh"

namespace Dune
{

  /**
   * \brief Layout map for Brezzi-Douglas-Marini-2 elements on quadrilaterals
   *
   * \ingroup LocalLayoutImplementation
   * \nosubgrouping
   * \implements Dune::LocalCoefficientsVirtualImp
   */
  class BDM2Cube2DLocalCoefficients
  {

  public:
    //! \brief Standard constructor
    BDM2Cube2DLocalCoefficients() : li(14)
    {
      for (std::size_t i = 0; i < 4; ++i)
      {
        li[3 * i] = LocalKey(i,1,0);
        li[3 * i + 1] = LocalKey(i,1,1);
        li[3 * i + 2] = LocalKey(i,1,2);
      }
      li[12] = LocalKey(0,0,0);
      li[13] = LocalKey(0,0,1);
    }

    //! \brief number of coefficients
    std::size_t size() const
    {
      return 14;
    }

    //! \brief get i'th index
    const LocalKey& localKey(std::size_t i) const
    {
      return li[i];
    }

  private:
    std::vector<LocalKey> li;
  };
} // end namespace Dune
#endif // DUNE_LOCALFUNCTIONS_BREZZIDOUGLASMARINI2_CUBE2D_LOCALCOEFFICIENTS_HH
