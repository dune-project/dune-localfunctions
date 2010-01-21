// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P2_3DLOCALCOEFFICIENTS_HH
#define DUNE_P2_3DLOCALCOEFFICIENTS_HH

#include <cstddef>
#include <vector>

#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{

  /**@ingroup LocalLayoutImplementation
         \brief Layout map for P23D elements

         \nosubgrouping
     \implements Dune::LocalCoefficientsVirtualImp
   */
  class P23DLocalCoefficients
  {
  public:
    //! \brief Standard constructor
    P23DLocalCoefficients () : li(10)
    {
      // Vertex shape functions
      li[0] = LocalKey(0,3,0);
      li[1] = LocalKey(1,3,0);
      li[2] = LocalKey(2,3,0);
      li[3] = LocalKey(3,3,0);

      // Edge bubbles
      li[4] = LocalKey(0,2,0);
      li[5] = LocalKey(2,2,0);
      li[6] = LocalKey(1,2,0);
      li[7] = LocalKey(3,2,0);
      li[8] = LocalKey(4,2,0);
      li[9] = LocalKey(5,2,0);
    }

    //! number of coefficients
    std::size_t size () const
    {
      return 10;
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
