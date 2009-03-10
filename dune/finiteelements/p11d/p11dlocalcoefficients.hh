// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_P11DLOCALCOEFFICIENTS_HH
#define DUNE_P11DLOCALCOEFFICIENTS_HH

#include <iostream>
#include <vector>

#include "../common/localcoefficients.hh"

namespace Dune
{

  /**@ingroup LocalLayoutImplementation
         \brief Layout map for P0 elements

         \nosubgrouping
   */
  class P11DLocalCoefficients : public LocalCoefficientsInterface<P11DLocalCoefficients>
  {
  public:
    //! \brief Standard constructor
    P11DLocalCoefficients () : li(2)
    {
      for (int i=0; i<2; i++)
        li[i] = LocalKey(i,1,0);
    }

    //! number of coefficients
    int size () const
    {
      return 2;
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

#endif
