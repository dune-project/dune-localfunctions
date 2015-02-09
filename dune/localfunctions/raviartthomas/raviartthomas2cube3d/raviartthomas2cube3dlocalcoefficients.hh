// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS2_CUBE3D_LOCALCOEFFICIENTS_HH
#define DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS2_CUBE3D_LOCALCOEFFICIENTS_HH

#include <cstddef>
#include <vector>

#include "../../common/localkey.hh"

namespace Dune
{

  /**
   * \ingroup LocalLayoutImplementation
   * \brief Layout map for Raviart-Thomas-2 elements on hexahedron.
   *
   * \nosubgrouping
   * \implements Dune::LocalCoefficientsVirtualImp
   */
  class RT2Cube3DLocalCoefficients
  {

  public:
    //! \brief Standard constructor
    RT2Cube3DLocalCoefficients () : li(108)
    {
      for (std::size_t i = 0; i < 6; i++)
      {
        li[9*i] = LocalKey(i,1,0);
        li[9*i + 1] = LocalKey(i,1,1);
        li[9*i + 2] = LocalKey(i,1,2);
        li[9*i + 3] = LocalKey(i,1,3);
        li[9*i + 4] = LocalKey(i,1,4);
        li[9*i + 5] = LocalKey(i,1,5);
        li[9*i + 6] = LocalKey(i,1,6);
        li[9*i + 7] = LocalKey(i,1,7);
        li[9*i + 8] = LocalKey(i,1,8);
      }

      for (std::size_t i = 0; i < 54; i++)
      {
        li[i + 54] = LocalKey(0,0,i);
      }
    }

    //! \brief number of coefficients
    std::size_t size () const
    {
      return 108;
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
#endif // DUNE_LOCALFUNCTIONS_RAVIARTTHOMAS2_CUBE3D_LOCALCOEFFICIENTS_HH
