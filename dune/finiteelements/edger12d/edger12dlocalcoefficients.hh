// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_EDGER12DLOCALCOEFFICIENTS_HH
#define DUNE_EDGER12DLOCALCOEFFICIENTS_HH

#include <iostream>
#include <vector>

#include "../common/localcoefficients.hh"

namespace Dune
{

  /**@ingroup LocalLayoutImplementation
     \brief Layout map for edge R1 elements

     \nosubgrouping
   */
  class EdgeR12DLocalCoefficients
    : public LocalCoefficientsInterface<EdgeR12DLocalCoefficients>
  {
  public:
    //! \brief Standard constructor
    EdgeR12DLocalCoefficients () : li(4)
    {
      li[0] = LocalKey(2,1,0);
      li[1] = LocalKey(3,1,0);
      li[2] = LocalKey(0,1,0);
      li[3] = LocalKey(1,1,0);
    }

    //! number of coefficients
    int size () const
    {
      return 4;
    }

    //! get i'th index
    const LocalKey& localKey (int i) const
    {
      return li[i];
    }

  private:
    std::vector<LocalKey> li;
  };

}

#endif // DUNE_EDGER12DLOCALCOEFFICIENTS_HH
