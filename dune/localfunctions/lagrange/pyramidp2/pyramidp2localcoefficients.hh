// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PYRAMID_P2_LOCALCOEFFICIENTS_HH
#define DUNE_PYRAMID_P2_LOCALCOEFFICIENTS_HH

#include <cstddef>
#include <vector>

#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{

  /**@ingroup LocalLayoutImplementation
     \brief Layout map for PyramidP2 elements

     \nosubgrouping
     \implements Dune::LocalCoefficientsVirtualImp
   */
  class PyramidP2LocalCoefficients
  {
  public:
    //! \brief Standard constructor
    PyramidP2LocalCoefficients () : li(14)
    {
      // Vertex shape functions
      li[0] = LocalKey(0,3,0);
      li[1] = LocalKey(1,3,0);
      li[2] = LocalKey(2,3,0);
      li[3] = LocalKey(3,3,0);
      li[4] = LocalKey(4,3,0);

      // Edge shape functions
      li[5] = LocalKey(0,2,0);
      li[6] = LocalKey(1,2,0);
      li[7] = LocalKey(2,2,0);
      li[8] = LocalKey(3,2,0);
      li[9] = LocalKey(4,2,0);
      li[10] = LocalKey(5,2,0);
      li[11] = LocalKey(6,2,0);
      li[12] = LocalKey(7,2,0);

      // base face shape function
      li[13] = LocalKey(0,1,0);
    }

    //! number of coefficients
    std::size_t size () const
    {
      return 14;
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
