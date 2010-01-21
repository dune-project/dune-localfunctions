// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PRISM_P2_LOCALCOEFFICIENTS_HH
#define DUNE_PRISM_P2_LOCALCOEFFICIENTS_HH

#include <cstddef>
#include <vector>

#include <dune/localfunctions/common/localkey.hh>

namespace Dune
{

  /**@ingroup LocalLayoutImplementation
         \brief Layout map for PrismP2 elements

         \nosubgrouping
     \implements Dune::LocalCoefficientsVirtualImp
   */
  class PrismP2LocalCoefficients
  {
  public:
    //! \brief Standard constructor
    PrismP2LocalCoefficients () : li(18)
    {
      // Vertex shape functions
      li[0] = LocalKey(0,3,0);
      li[1] = LocalKey(1,3,0);
      li[2] = LocalKey(2,3,0);
      li[3] = LocalKey(3,3,0);
      li[4] = LocalKey(4,3,0);
      li[5] = LocalKey(5,3,0);

      // Edge shape functions
      li[6] = LocalKey(0,2,0);
      li[7] = LocalKey(1,2,0);
      li[8] = LocalKey(2,2,0);
      li[9] = LocalKey(3,2,0);
      li[10] = LocalKey(4,2,0);
      li[11] = LocalKey(5,2,0);
      li[12] = LocalKey(6,2,0);
      li[13] = LocalKey(7,2,0);
      li[14] = LocalKey(8,2,0);

      // Quadrilateral sides shape functions
      li[15] = LocalKey(0,1,0);
      li[16] = LocalKey(1,1,0);
      li[17] = LocalKey(2,1,0);
    }

    //! number of coefficients
    std::size_t size () const
    {
      return 18;
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
